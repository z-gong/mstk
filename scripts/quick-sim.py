#!/usr/bin/env python3

import sys
import pickle
import traceback
import argparse
import numpy as np
from openmm import openmm as mm
from mstk.topology import Topology, Molecule
from mstk.trajectory import Trajectory
from mstk.forcefield import ForceField
from mstk.forcefield.typer import ZftTyper, typer_primitive
from mstk.forcefield.errors import *
from mstk.simsys import System
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
    parser.add_argument('--box', nargs='+', type=float,
                        help='periodic box size if not provided by topology or configuration')
    parser.add_argument('--packmol', action='store_true',
                        help='generate Packmol input files for building coordinates')
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
        top = Topology.open(inp)
        ref_mols += top.molecules

    if args.number is None:
        args.number = [1] * len(ref_mols)

    if len(args.number) != len(ref_mols):
        raise Exception('inconsistent molecules and numbers')

    molecules = []
    for mol, n in zip(ref_mols, args.number):
        molecules += [mol] * n
    top = Topology(molecules)
    top.remove_drude_particles()
    top.remove_virtual_sites()

    if args.conf is not None:
        frame = Trajectory.read_frame_from_file(args.conf, -1)
        if frame.cell.volume != 0:
            top.cell.set_box(frame.cell.vectors)

    if args.box is not None:
        if len(args.box) == 3:
            box = args.box
        elif len(args.box) == 1:
            box = args.box * 3
        else:
            raise Exception('box should be float or list of three floats')
        top.cell.set_box(box)

    ff = ForceField.open(*args.ff)

    mol_count = dict(zip(ref_mols, args.number))
    for mol in mol_count.keys():
        logger.info('Processing %s ...' % str(mol))
        if typer is not None:
            try:
                typer.type(mol)
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
                mol.guess_connectivity_from_ff(ff, pbc='xyz', cell=top.cell)
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
