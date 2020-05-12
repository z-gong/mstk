#!/usr/bin/env python3

import argparse
import warnings
from mstools.topology import Topology
from mstools.trajectory import Trajectory
from mstools.forcefield import FFSet, ZftTyper
from mstools.simsys import System

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', nargs='+', required=True, type=str,
                    help='Topology files for molecules. '
                         'String starts with : will be treated as SMILES')
parser.add_argument('-f', '--forcefield', nargs='+', required=True, type=str,
                    help='Forcefield files')
parser.add_argument('-n', '--number', nargs='+', type=int, help='number of molecules')
parser.add_argument('--typer', type=str,
                    help='Typing file. Required if SMILES provided for topology')
parser.add_argument('--trj', type=str,
                    help='Trajectory file for positions and box')
parser.add_argument('--box', nargs=3, type=float,
                    help='Periodic box size if not provided by topology or trajectory')
parser.add_argument('--packmol', action='store_true',
                    help='Generate Packmol input files for building coordinates')
args = parser.parse_args()

if args.number is None:
    args.number = [1] * len(args.input)

if len(args.input) != len(args.number):
    raise Exception('input and number inconsistent')

if args.typer is not None:
    typer = ZftTyper(args.typer)
else:
    typer = None

molecules = []
for inp, n in zip(args.input, args.number):
    top = Topology.open(inp)
    if inp.startswith(':'):
        if typer is None:
            raise Exception('--typer is required for SMILES input')
        typer.type(top)
    molecules += top.molecules * n
top = Topology(molecules)

if args.trj is not None:
    trj = Trajectory.open(args.trj)
    frame = trj.read_frame(trj.n_frame - 1)
    top.set_positions(frame.positions)
    if frame.cell.volume != 0:
        top.cell.set_box(frame.cell.vectors)

if args.box is not None:
    top.cell.set_box(args.box)

ff = FFSet.open(*args.forcefield)
if len(ff.polarizable_terms) > 0:
    top.generate_drude_particles(ff)
top.assign_mass_from_ff(ff)
top.assign_charge_from_ff(ff)

system = System(top, ff)
if args.packmol:
    mol_numbers = top.get_unique_molecules()
    top.update_molecules(list(mol_numbers.keys()))
    top.scale_with_packmol(list(mol_numbers.values()))
else:
    if args.trj is None:
        warnings.warn('Trajectory file not provided, '
                      'will use the positions and cell from the topology')
    try:
        system.export_gmx(gro_out='_conf.gro', top_out='_topol.top', mdp_out='_grompp.mdp')
    except Exception as e:
        print(repr(e))
    top.write('_topol.psf')
