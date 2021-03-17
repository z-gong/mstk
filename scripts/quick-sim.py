#!/usr/bin/env python3

import sys
import argparse
import numpy as np
from mstools.topology import Topology, UnitCell, Molecule
from mstools.trajectory import Trajectory
from mstools.forcefield import ForceField, ZftTyper, PaduaLJScaler
from mstools.forcefield.errors import *
from mstools.simsys import System
from mstools import logger

parser = argparse.ArgumentParser()
parser.add_argument('input', nargs='+', type=str,
                    help='topology files for molecules. '
                         'String starts with : will be treated as SMILES')
parser.add_argument('-f', '--forcefield', nargs='+', required=True, type=str,
                    help='forcefield files')
parser.add_argument('-n', '--number', nargs='+', type=int, help='number of molecules')
parser.add_argument('--typer', type=str,
                    help='typing file. Required if SMILES provided for topology')
parser.add_argument('--ljscale', type=str, help='input file for empirical LJ scaling')
parser.add_argument('--nodrude', action='store_true',
                    help='disable Drude particles and virtual sites even if present in force field')
parser.add_argument('--trj', type=str,
                    help='trajectory file for positions and box. The last frame will be used')
parser.add_argument('--box', nargs=3, type=float,
                    help='periodic box size if not provided by topology or trajectory')
parser.add_argument('--packmol', action='store_true',
                    help='generate Packmol input files for building coordinates')
parser.add_argument('--qscale', default=1, type=float, help='scale the charge of atoms')
parser.add_argument('--scaleeps', type=float, default=1.0,
                    help='extra scaling parameter for all LJ epsilon')
parser.add_argument('--scalesig', type=float, default=1.0,
                    help='extra scaling parameter for all LJ sigma')
parser.add_argument('--scaleignoreatom', nargs='+', type=str, default=[],
                    help='ignore these atom types for extra scaling')
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
    if inp.startswith(':'):
        if typer is None:
            raise Exception('--typer is required for SMILES input')
    top = Topology.open(inp)
    molecules += top.molecules * n
top = Topology(molecules)
top.remove_drude_particles()
top.remove_virtual_sites()

if args.trj is not None:
    frame = Trajectory.read_frame_from_file(args.trj, -1)
    if frame.cell.volume != 0:
        top.cell.set_box(frame.cell.vectors)

if args.box is not None:
    top.cell.set_box(args.box)

ff = ForceField.open(*args.forcefield)
if args.nodrude:
    ff.polarizable_terms.clear()
    ff.virtual_site_terms.clear()
if args.ljscale is not None:
    scaler = PaduaLJScaler(args.ljscale)
    scaler.scale(ff)

if args.scaleeps != 1.0 or args.scalesig != 1.0:
    for vdw in list(ff.vdw_terms.values()) + list(ff.pairwise_vdw_terms.values()):
        if vdw.type1 in args.scaleignoreatom or vdw.type2 in args.scaleignoreatom:
            continue
        if args.scaleeps != 1.0:
            vdw.epsilon *= args.scaleeps
            vdw.comments.append('eps*%.3f' % args.scaleeps)
        if args.scalesig != 1.0:
            vdw.sigma *= args.scalesig
            vdw.comments.append('sig*%.3f' % args.scalesig)

mol_count = top.get_unique_molecules(deepcopy=False)
for mol in mol_count.keys():
    mol: Molecule
    logger.info('Processing %s ...' % str(mol))
    if typer is not None:
        try:
            typer.type_molecule(mol)
        except TypingNotSupportedError as e:
            pass
        except TypingUndefinedError as e:
            xyz = '_typing_' + mol.name + '.xyz'
            Topology([mol]).write(xyz)
            logger.error('Failed typing %s: %s. Check %s' % (mol, str(e), xyz))
            sys.exit(1)

    if mol.n_atom > 1 and mol.n_bond == 0:
        if mol.has_position:
            logger.warning(f'{str(mol)} carries no bond. Guessing connectivity from FF')
            mol.guess_connectivity_from_ff(ff, angle_tolerance=15, pbc='xyz', cell=top.cell)
        else:
            logger.warning(f'{str(mol)} carries no bond. Make sure the topology is correct')

    if ff.is_polarizable:
        mol.generate_drude_particles(ff)
    if ff.has_virtual_site:
        mol.generate_virtual_sites(ff)
    mol.assign_mass_from_ff(ff)
    mol.assign_charge_from_ff(ff)

if args.packmol:
    top.update_molecules(list(mol_count.keys()))
    top.scale_with_packmol(list(mol_count.values()))
    sys.exit(0)

logger.info('Exporting ...')

top.update_molecules(list(mol_count.keys()), list(mol_count.values()))

if args.trj is None:
    logger.warning('Trajectory file not provided. '
                   'Will use positions and cell from the topology')
else:
    _positions_set = False
    if len(frame.positions) == top.n_atom:
        top.set_positions(frame.positions)
        _positions_set = True
    elif ff.is_polarizable or ff.has_virtual_site:
        atoms_real = [atom for atom in top.atoms if not atom.is_drude and atom.virtual_site is None]
        if len(frame.positions) == len(atoms_real):
            for i, atom in enumerate(atoms_real):
                atom.position = frame.positions[i][:]
            np.random.seed(1)
            for parent, drude in top.get_drude_pairs():
                drude.position = parent.position + (np.random.random(3) - 0.5) / 100
            for parent, vsite in top.get_virtual_site_pairs():
                vsite.position = vsite.virtual_site.calc_position()
            _positions_set = True

    if not _positions_set:
        logger.error(
            'Numbers of atoms in trajectory (%i) and topology (all: %i, real: %i) do not match' % (
                len(frame.positions), top.n_atom,
                top.n_atom - len(top.get_drude_pairs()) - len(top.get_virtual_site_pairs())))
        sys.exit(1)

for atom in top.atoms:
    if args.qscale != 1 and atom.type not in args.scaleignoreatom:
        atom.charge *= args.qscale

system = System(top, ff)
system.export_gromacs(gro_out='_conf.gro', top_out='_topol.top', mdp_out='_grompp.mdp')
system.export_charmm(pdb_out=None, psf_out='_topol.psf', prm_out='_ff.prm')
system.export_lammps(data_out='_data.lmp', in_out='_in.lmp')
