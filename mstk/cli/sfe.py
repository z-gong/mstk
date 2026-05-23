#!/usr/bin/env python3

import os
import argparse
import numpy as np
from openmm import openmm as mm, app
from mstk import logger
from mstk.ommhelper import get_platform_properties, unit
from mstk.chem import constant
from mstk.topology import Topology
from mstk.trajectory import Trajectory
from mstk.forcefield import ForceField
from mstk.simsys import System
from mstk.sfe import SFEManager


def add_subcommand(subparsers):
    parser = subparsers.add_parser('sfe', help='Solvation free energy via alchemical decoupling',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub = parser.add_subparsers(dest='sfe_command', required=True)

    # sfe run
    run_parser = sub.add_parser('run', help='Run a single SFE lambda window',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                allow_abbrev=False)
    run_parser.add_argument('-p', '--top', type=str, required=True,
                            help='topology file')
    run_parser.add_argument('-c', '--conf', type=str, required=True,
                            help='configuration file')
    run_parser.add_argument('-f', '--ff', nargs='+', type=str, required=True,
                            help='force field file(s)')
    run_parser.add_argument('--molid', type=int, required=True,
                            help='index of the molecule to decouple')
    run_parser.add_argument('-w', '--window', type=int, required=True,
                            help='lambda window index (0-based)')
    run_parser.add_argument('--nwindow', type=int, default=16,
                            help='total number of lambda windows')
    run_parser.add_argument('-n', '--nstep', type=int, default=500000,
                            help='production steps')
    run_parser.add_argument('--nstepeq', type=int, default=100000,
                            help='equilibration steps')
    run_parser.add_argument('-t', '--temp', type=float, default=300.0,
                            help='temperature in K')
    run_parser.add_argument('--press', type=float, default=1.0,
                            help='pressure in bar')
    run_parser.add_argument('--dt', type=float, default=0.002,
                            help='timestep in ps')
    run_parser.add_argument('-o', '--output', type=str, default=None,
                            help='output CSV file (default: dU_<window>.csv)')
    run_parser.set_defaults(func=_run)

    # sfe mbar
    mbar_parser = sub.add_parser('mbar', help='Run MBAR analysis on SFE output files',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    mbar_parser.add_argument('input_files', nargs='+', type=str,
                             help='CSV files from each lambda window')
    mbar_parser.set_defaults(func=_mbar)

    parser.set_defaults(func=lambda args: parser.print_help())


def _run(args):

    output = args.output or f'dU_{args.window}.csv'

    for arg, val in vars(args).items():
        if arg not in ('func', 'command', 'sfe_command'):
            logger.info(f'--{arg:15s} {val}')

    top = Topology.open(args.top)
    frame = Trajectory.read_frame_from_file(args.conf, -1)
    top.cell.set_box(frame.cell.vectors)
    top.set_positions(frame.positions)

    ff = ForceField.open(*args.ff)
    system = System(top, ff)

    # Build alchemical system.
    sfe = SFEManager(system, args.molid, n_lambda=args.nwindow)
    sfe.omm_system.addForce(mm.MonteCarloBarostat(args.press * unit.bar, args.temp * unit.kelvin, 100))

    lambda_coul, lambda_vdw = sfe.schedule[args.window]
    logger.info(f'Window {args.window}/{args.nwindow}: '
                f'lambda_coul={lambda_coul:.4f}, lambda_vdw={lambda_vdw:.4f}')

    # Create simulation.
    platform, properties = get_platform_properties()
    integrator = mm.LangevinMiddleIntegrator(args.temp * unit.kelvin, 1.0 / unit.ps, args.dt * unit.ps)
    sim = app.Simulation(top.to_omm_topology(), sfe.omm_system, integrator, platform, properties)
    sim.context.setPositions(top.positions)

    # Set initial lambda state.
    sfe.set_lambda_state(sim.context, args.window)

    # Minimize.
    logger.info('Minimizing...')
    sim.minimizeEnergy(maxIterations=1000)

    # Equilibration.
    logger.info(f'Equilibrating for {args.nstepeq} steps...')
    sim.step(args.nstepeq)

    # Production: collect dU every 200 steps.
    sample_interval = 200
    n_samples = args.nstep // sample_interval
    n_lambda = args.nwindow
    logger.info(f'Production: {args.nstep} steps, collecting {n_samples} samples...')

    fout = open(output, 'w')
    fout.write(f'# window={args.window} nwindow={n_lambda} temp={args.temp} press={args.press}\n')
    fout.write(','.join([f'dU_{k}' for k in range(n_lambda)]) + '\n')

    for sample in range(n_samples):
        sim.step(sample_interval)

        U_current = sim.context.getState(getEnergy=True).getPotentialEnergy()._value

        dU_row = np.zeros(n_lambda)
        for k in range(n_lambda):
            if k == args.window:
                continue
            sfe.set_lambda_state(sim.context, k)
            U_k = sim.context.getState(getEnergy=True).getPotentialEnergy()._value
            dU_row[k] = U_k - U_current

        sfe.set_lambda_state(sim.context, args.window)
        fout.write(','.join(f'{v:.6e}' for v in dU_row) + '\n')
        fout.flush()

    fout.close()
    logger.info(f'Wrote {n_samples} samples to {output}')


def _parse_dU_header(filepath):
    with open(filepath) as f:
        line = f.readline()
    if not line.startswith('#'):
        raise ValueError(f'{filepath}: missing metadata header')
    meta = {}
    for token in line[1:].split():
        key, val = token.split('=')
        meta[key] = val
    for key in ('window', 'nwindow', 'temp', 'press'):
        if key not in meta:
            raise ValueError(f'{filepath}: missing "{key}" in header')
    return {
        'window': int(meta['window']),
        'nwindow': int(meta['nwindow']),
        'temp': float(meta['temp']),
        'press': float(meta['press']),
        'file': filepath,
    }


def _mbar(args):
    headers = [_parse_dU_header(f) for f in args.input_files]

    nwindows = set(h['nwindow'] for h in headers)
    temps = set(h['temp'] for h in headers)
    presses = set(h['press'] for h in headers)
    if len(nwindows) > 1:
        raise ValueError(f'Inconsistent nwindow across files: {nwindows}')
    if len(temps) > 1:
        raise ValueError(f'Inconsistent temp across files: {temps}')
    if len(presses) > 1:
        raise ValueError(f'Inconsistent press across files: {presses}')

    n_lambda = nwindows.pop()
    temp = temps.pop()

    windows = [h['window'] for h in headers]
    if len(set(windows)) != len(windows):
        raise ValueError(f'Duplicate window indices found')
    if set(windows) != set(range(n_lambda)):
        missing = set(range(n_lambda)) - set(windows)
        raise ValueError(f'Missing windows: {sorted(missing)}')

    headers.sort(key=lambda h: h['window'])
    logger.info(f'Found {n_lambda} SFE windows, temp={temp} K')

    R = constant.BOLTZMANN * constant.AVOGADRO / 1000  # kJ/(mol*K)
    beta = 1.0 / (R * temp)

    all_data = []
    for h in headers:
        data = np.loadtxt(h['file'], delimiter=',', skiprows=2)
        all_data.append(data)

    N_k = np.array([d.shape[0] for d in all_data])
    n_total = int(N_k.sum())

    u_kn = np.zeros((n_lambda, n_total))

    offset = 0
    for data in all_data:
        n_i = data.shape[0]
        for k in range(n_lambda):
            u_kn[k, offset:offset + n_i] = beta * data[:, k]
        offset += n_i

    mu_ex, error, dG_adj, dG_adj_err = SFEManager.compute_mbar(u_kn, N_k, temp)

    logger.info('Per-window dG (kJ/mol):')
    for i in range(len(dG_adj)):
        logger.info(f'  {i:2d} -> {i + 1:2d}: {dG_adj[i]:8.3f} +/- {dG_adj_err[i]:.3f}')
    logger.info(f'MBAR result: mu_ex = {mu_ex:.3f} +/- {error:.3f} kJ/mol')
