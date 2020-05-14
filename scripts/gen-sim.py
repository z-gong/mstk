#!/usr/bin/env python3

import sys
import argparse
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
parser.add_argument('--ljscale', type=str, help='input files for empirical LJ scaling')
parser.add_argument('--nodrude', action='store_true',
                    help='disable Drude particles even if it is polarizable force field')
parser.add_argument('--trj', type=str,
                    help='trajectory file for positions and box. The last frame will be used')
parser.add_argument('--box', nargs=3, type=float,
                    help='periodic box size if not provided by topology or trajectory')
parser.add_argument('--packmol', action='store_true',
                    help='generate Packmol input files for building coordinates')
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
if args.nodrude:
    top.remove_drude_particles()
if args.box is not None:
    top.cell.set_box(args.box)

ff = ForceField.open(*args.forcefield)
if args.nodrude:
    ff.polarizable_terms = {}
if args.ljscale is not None:
    scaler = PaduaLJScaler(args.ljscale)
    scaler.scale(ff)
    logger.warning('LJ scaling file provided. Check the generated FF carefully')

unique_molecules = top.get_unique_molecules()
for mol in unique_molecules.keys():
    mol: Molecule
    logger.info('Processing %s ...' % str(mol))
    try:
        typer.type_molecule(mol)
    except TypingNotSupportedError as e:
        pass
    except TypingUndefinedError as e:
        logger.error('Typing failed %s: %s' % (mol, str(e)))

    if mol.n_bond == 0 and mol.has_position:
        logger.warning(f'{str(mol)} carries no bond. Try to guess connectivity from forcefield')
        mol.guess_connectivity_from_ff(ff, angle_tolerance=15, pbc='xyz', cell=top.cell)

    if ff.is_polarizable:
        mol.generate_drude_particles(ff)
    mol.assign_mass_from_ff(ff)
    mol.assign_charge_from_ff(ff)

top.update_molecules(list(unique_molecules.keys()), list(unique_molecules.values()))

if args.trj is not None:
    frame = Trajectory.read_frame_from_file(args.trj, -1)
    if len(frame.positions) != top.n_atom:
        logger.error('Number of atoms in trjectory and topology doesn\'t match')
        sys.exit(1)

    top.set_positions(frame.positions)
    if args.box is None and frame.cell.volume != 0:
        top.cell.set_box(frame.cell.vectors)

logger.info('Exporting ...')
if args.packmol:
    mol_numbers = top.get_unique_molecules()
    top.update_molecules(list(mol_numbers.keys()))
    top.scale_with_packmol(list(mol_numbers.values()))
else:
    if args.trj is None:
        logger.warning('Trajectory file not provided. '
                       'Will use positions and cell from the topology')
    system = System(top, ff)
    system.export_gromacs(gro_out='_conf.gro', top_out='_topol.top', mdp_out='_grompp.mdp')
    system.export_charmm(pdb_out=None, psf_out='_topol.psf', prm_out='_ff.prm')
