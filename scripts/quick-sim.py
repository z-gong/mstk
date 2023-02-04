#!/usr/bin/env python3

import sys
import pickle
import traceback
import argparse
import random
import numpy as np
from openmm import openmm as mm
from mstk.topology import Topology, Molecule
from mstk.trajectory import Trajectory
from mstk.forcefield import ForceField
from mstk.forcefield.typer import ZftTyper, typer_primitive
from mstk.forcefield.errors import *
from mstk.simsys import System
from mstk.chem import constant
from mstk import logger


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs='+', type=str,
                        help='topology files for molecules. '
                             'String starts with : will be treated as SMILES')
    parser.add_argument('-n', '--number', nargs='+', type=int, help='number of molecules')
    parser.add_argument('-f', '--ff', nargs='+', required=True, type=str,
                        help='forcefield files')
    parser.add_argument('-t', '--typer', type=str,
                        help='typing file. Required if SMILES provided for topology. '
                             'If set to primitive, then the ZftTyper defined by primitive.zft will be used')
    parser.add_argument('-c', '--conf', type=str,
                        help='configuration file with positions and box. The last frame will be used')
    parser.add_argument('--density', type=float,
                        help='the density of the target box in g/mL. '
                             'The size of cubic box will be calculated from it. '
                             'If the argument box is also provided, density will be ignored')
    parser.add_argument('--box', nargs='+', type=float,
                        help='periodic box size if not provided by topology or configuration')
    parser.add_argument('--packmol', action='store_true',
                        help='generate Packmol input files for building coordinates')
    parser.add_argument('--ua', action='store_true', help='export united-atom model')
    return parser.parse_args()


def main():
    args = parse_args()
    for arg, val in vars(args).items():
        logger.info(f'--{arg:10s} {val}')

    if args.typer is not None:
        if args.typer == 'primitive':
            typer = typer_primitive
        else:
            typer = ZftTyper(args.typer)
    else:
        typer = None

    ref_mols = []
    for inp in args.input:
        if inp.startswith(':') or inp.lower().endswith('.smi'):
            if typer is None:
                raise Exception('--typer is required for SMILES input')
        ref_mols += Topology.open(inp).molecules

    if args.number is None:
        args.number = [1] * len(ref_mols)

    if len(args.number) != len(ref_mols):
        raise Exception('Inconsistent molecules and numbers')

    # in case number is 0 for some molecules
    numbers = list(args.number)
    for i in range(len(numbers) - 1, -1, -1):
        if args.number[i] <= 0:
            ref_mols.pop(i)
            numbers.pop(i)

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

    top = Topology(ref_mols, numbers)

    if args.conf is not None:
        frame = Trajectory.read_frame_from_file(args.conf, -1)
        # do not set position because the frame may not contain Drude particles and virtual sites
        if frame.cell.volume != 0:
            top.cell.set_box(frame.cell.vectors)

    if args.density is not None:
        mass = sum(atom.mass for atom in top.atoms)  # g/mol
        vol = mass / args.density / constant.AVOGADRO  # cm^3
        box = [vol ** (1 / 3) * constant.CENTI / constant.NANO] * 3  # nm
        top.cell.set_box(box)

    if args.box is not None:
        if len(args.box) == 3:
            box = args.box
        elif len(args.box) == 1:
            box = args.box * 3
        else:
            raise Exception('box should be float or list of three floats')
        top.cell.set_box(box)

    if args.packmol:
        top.update_molecules(ref_mols)
        top.scale_with_packmol(numbers, seed=random.randrange(int(1E5), int(1E6)), tempdir='.')
        sys.exit(0)

    logger.info('Exporting ...')
    if args.conf is None:
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

    system = System(top, ff)
    logger.info('Exporting NAMD...')
    try:
        system.export_namd(pdb_out='_namd-conf.pdb', psf_out='_namd-topol.psf', prm_out=None, pdb_atom_type=True)
    except Exception as e:
        logger.error('Failed exporting NAMD')
        traceback.print_exc()

    logger.info('Exporting GROMACS...')
    try:
        system.export_gromacs(gro_out='_gmx-conf.gro', top_out='_gmx-topol.top', mdp_out='_gmx-pp.mdp')
    except Exception as e:
        logger.error('Failed exporting GROMACS')
        traceback.print_exc()

    logger.info('Exporting LAMMPS...')
    try:
        system.export_lammps(data_out='_lmp-data.lmp', in_out='_lmp-in.lmp')
    except Exception as e:
        logger.error('Failed exporting LAMMPS')
        traceback.print_exc()

    logger.info('Exporting OpenMM...')
    try:
        omm_top = top.to_omm_topology()
        with open('_omm-top.pkl', 'wb') as f:
            f.write(pickle.dumps(omm_top))
        omm_sys = system.to_omm_system()
        with open('_omm-sys.xml', 'w') as f:
            f.write(mm.XmlSerializer.serialize(omm_sys))
    except Exception as e:
        logger.error('Failed exporting OpenMM')
        traceback.print_exc()


if __name__ == '__main__':
    main()
