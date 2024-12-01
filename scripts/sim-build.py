#!/usr/bin/env python3

import sys
import argparse
import random
import numpy as np
from mstk.topology import Topology, Molecule
from mstk.trajectory import Trajectory, Frame
from mstk.forcefield import ForceField
from mstk.forcefield.typer import Typer, GaffTyper
from mstk.forcefield.errors import *
from mstk.chem import constant
from mstk import logger


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--top', nargs='+', type=str, required=True,
                        help='topology files for molecules. '
                             'String starts with : will be treated as SMILES')
    parser.add_argument('-c', '--conf', nargs='+', type=str,
                        help='configuration file for each molecule. The last frame will be used. '
                             'Positions are required to build box and write coordinates. '
                             'If the topology doesn\'t contain positions, it must be provided by this argument. '
                             'The number of conf files can be smaller or equal to the number of molecule types')
    parser.add_argument('--op', type=str, default='top.psf', help='output topology file')
    parser.add_argument('--oc', type=str, default='conf.gro', help='output configuration file')
    parser.add_argument('-n', '--number', nargs='+', type=int, help='number of each molecule')
    parser.add_argument('-f', '--ff', nargs='+', required=True, type=str, help='forcefield files')
    parser.add_argument('-t', '--typer', type=str,
                        help='typing file. Required if SMILES provided for topology. '
                             'If set to gaff, will use the GAFF typing engine')
    parser.add_argument('--ua', action='store_true',
                        help='use united-atom model. Remove non-polar hydrogens (connected to C/Si) from molecules')
    parser.add_argument('--box', nargs='+', type=float,
                        help='periodic box size. '
                             'The box will be rectangular if three values are provided. '
                             'The box will be cubic if only one value is provided.'
                        )
    parser.add_argument('--density', type=float,
                        help='the density of the target box in g/mL. '
                             'The size of cubic box will be calculated from it. '
                             'If the argument box is provided, density will be ignored')
    parser.add_argument('--packmol', action='store_true',
                        help='generate Packmol input files instead of write coordinates directly')
    return parser.parse_args()


def main():
    args = parse_args()
    for arg, val in vars(args).items():
        logger.info(f'--{arg:10s} {val}')

    if args.typer == 'gaff':
        typer = GaffTyper()
    elif args.typer:
        typer = Typer.open(args.typer)
    else:
        typer = None

    ref_mols = []
    for inp in args.top:
        if inp.startswith(':') or inp.lower().endswith('.smi'):
            if typer is None:
                raise Exception('--typer is required for SMILES input')
        ref_mols += Topology.open(inp).molecules

    args.number = args.number or [1] * len(ref_mols)
    if len(args.number) != len(ref_mols):
        raise Exception('Inconsistent molecules and numbers')

    ff = ForceField.open(*args.ff)

    for mol in ref_mols:
        logger.info(f'Processing {mol}...')
        mol.remove_drude_particles()
        mol.remove_virtual_sites()
        if mol.n_atom > 1 and mol.n_bond == 0:
            logger.warning(f'{mol} carries no bond. Make sure the topology is correct')

        if typer is not None:
            try:
                typer.type(mol)
            except TypingNotSupportedError as e:
                logger.warning(f'Failed typing {mol}: {e}. '
                               f'This warning can be ignored if atom type is already assigned')
            except TypingUndefinedError as e:
                xyz = '_typing_' + mol.name + '.xyz'
                Topology([mol]).write(xyz)
                logger.error(f'Failed typing {mol}: {e}. Check {xyz}')
                sys.exit(1)

        if args.ua:
            mol.remove_non_polar_hydrogens()

        if ff.is_polarizable:
            mol.generate_drude_particles(ff)
        if ff.has_virtual_site:
            mol.generate_virtual_sites(ff)
        ff.assign_mass(mol)
        ff.assign_charge(mol)

    conf_files = args.conf or []
    for mol, conf_file in zip(ref_mols, conf_files):
        frame = Trajectory.read_frame_from_file(conf_file, -1)
        _positions_set = False
        if len(frame.positions) == mol.n_atom:
            mol.set_positions(frame.positions)
            _positions_set = True
        elif ff.is_polarizable or ff.has_virtual_site:
            atoms_real = [atom for atom in mol.atoms if not atom.is_drude and atom.virtual_site is None]
            if len(frame.positions) == len(atoms_real):
                for i, atom in enumerate(atoms_real):
                    atom.position = frame.positions[i][:]
                np.random.seed(1)
                for parent, drude in mol.get_drude_pairs():
                    drude.position = parent.position + (np.random.random(3) - 0.5) / 100
                for parent, vsite in mol.get_virtual_site_pairs():
                    vsite.position = vsite.virtual_site.calc_position()
                _positions_set = True

        if not _positions_set:
            logger.error(f'Numbers of atoms in {conf_file} (%i) and {mol} (all: %i, real: %i) do not match' % (
                len(frame.positions), mol.n_atom,
                mol.n_atom - len(mol.get_drude_pairs()) - len(mol.get_virtual_site_pairs())))
            sys.exit(1)

    top = Topology(ref_mols, args.number)

    box = [0.0, 0.0, 0.0]
    if args.density is not None:
        mass = sum(atom.mass for atom in top.atoms)  # g/mol
        vol = mass / args.density / constant.AVOGADRO  # cm^3
        box = [vol ** (1 / 3) * constant.CENTI / constant.NANO] * 3  # nm

    if args.box is not None:
        if len(args.box) == 3:
            box = args.box
        elif len(args.box) == 1:
            box = args.box * 3
        else:
            raise Exception('box should be one or three floats')

    logger.info(f'box = {box}')
    top.cell.set_box(box)

    top.write(args.op)

    if args.packmol:
        top.update_molecules(ref_mols)
        # write the Packmol input files only. Can modify it before execution
        top.scale_with_packmol(args.number, seed=random.randrange(int(1E5), int(1E6)), tempdir='.')
    else:
        frame = Frame(top.n_atom)
        frame.cell = top.cell
        frame.positions = top.positions
        trj = Trajectory.open(args.oc, 'w')
        trj.write_frame(frame, top)
        trj.close()


if __name__ == '__main__':
    main()
