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
from mstk.simsys import System
from mstk.chem import constant
from mstk import logger


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--top', type=str, required=True,
                        help='topology file for the system. '
                             'This topology must have charges and masses assigned. '
                             'The charge and mass from force field will not be used')
    parser.add_argument('-f', '--ff', nargs='+', type=str, required=True, help='forcefield files')
    parser.add_argument('-c', '--conf', type=str, required=True,
                        help='configuration file with positions and box. The last frame will be used')
    parser.add_argument('--box', nargs='+', type=float,
                        help='periodic box size. '
                             'The box will be rectangular if three values are provided. '
                             'The box will be cubic if only one value is provided.'
                        )
    parser.add_argument('--density', type=float,
                        help='the density of the target box in g/mL. '
                             'The size of cubic box will be calculated from it. '
                             'If the argument box is provided, density will be ignored')
    parser.add_argument('--lammps', action='store_true', help='export LAMMPS files')
    parser.add_argument('--gromacs', action='store_true', help='export GROMACS files')
    parser.add_argument('--namd', action='store_true', help='export NAMD files')
    parser.add_argument('--openmm', action='store_true', help='export OpenMM files')
    parser.add_argument('--all', action='store_true', help='export files for all supported MD engines')
    return parser.parse_args()


def main():
    args = parse_args()
    for arg, val in vars(args).items():
        logger.info(f'--{arg:10s} {val}')

    top = Topology.open(args.top)
    ff = ForceField.open(*args.ff)

    frame = Trajectory.read_frame_from_file(args.conf, -1)
    if frame.cell.volume != 0:
        top.cell.set_box(frame.cell.vectors)

    # the frame may not contain Drude particles and virtual sites
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
        logger.error('Numbers of atoms in trajectory (%i) and topology (all: %i, real: %i) do not match' % (
            len(frame.positions), top.n_atom,
            top.n_atom - len(top.get_drude_pairs()) - len(top.get_virtual_site_pairs())))
        sys.exit(1)

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

    system = System(top, ff)

    if args.lammps or args.all:
        logger.info('Exporting LAMMPS...')
        try:
            system.export_lammps(data_out='_lmp-data.lmp', in_out='_lmp-in.lmp')
        except Exception as e:
            logger.error('Failed exporting LAMMPS')
            traceback.print_exc()

    if args.gromacs or args.all:
        logger.info('Exporting GROMACS...')
        try:
            system.export_gromacs(gro_out='_gmx-conf.gro', top_out='_gmx-topol.top', mdp_out='_gmx-pp.mdp')
        except Exception as e:
            logger.error('Failed exporting GROMACS')
            traceback.print_exc()

    if args.namd or args.all:
        logger.info('Exporting NAMD...')
        try:
            system.export_namd(pdb_out='_namd-conf.pdb', psf_out='_namd-top.psf', prm_out=None)
        except Exception as e:
            logger.error('Failed exporting NAMD')
            traceback.print_exc()

    if args.openmm or args.all:
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
