#!/usr/bin/env python3

import argparse
import numpy as np
from mstk.topology import Topology
from mstk.trajectory import Trajectory
from mstk.forcefield import ForceField
from mstk import logger


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--top', nargs='+', type=str, required=True, help='topology files')
    parser.add_argument('-n', '--number', nargs='+', type=int, help='number of molecules')
    parser.add_argument('-c', '--conf', type=str,
                        help='configuration file with positions and box. The last frame will be used. '
                             'This is for writing topology in formats supporting positions, e.g. PDB, XYZ')
    parser.add_argument('-o', '--output', required=True, type=str, help='output topology file')
    parser.add_argument('--ignore', nargs='+', default=[], type=str, help='ignore these molecule types')
    parser.add_argument('-f', '--ff', type=str, help='reassign charge from FF')
    parser.add_argument('--qscale', default=1, type=float, help='scale the charge of atoms')
    parser.add_argument('--qscaleignore', nargs='+', default=[], type=str,
                        help='ignore these molecule names for charge scaling')
    parser.add_argument('--box', nargs='+', type=float,
                        help='overwrite the box dimensions')
    parser.add_argument('--shift', nargs=3, default=[0, 0, 0], type=float,
                        help='shift the positions of all atoms')
    parser.add_argument('--ua', action='store_true',
                        help='remove hydrogen atoms bonded to C/Si. '
                             'Also remove the relevant bonds/angles/dihedrals/impropers')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    top_list = [Topology.open(inp) for inp in args.top]
    if args.number is None:
        args.number = [1] * len(top_list)

    molecules = []
    for top, n in zip(top_list, args.number):
        molecules.extend(top.molecules * n)
    if args.ignore:
        molecules = [mol for mol in molecules if mol.name not in args.ignore]

    top = top_list[0]
    top.update_molecules(molecules)
    logger.info(str(top))

    if args.conf:
        frame = Trajectory.read_frame_from_file(args.conf, -1)
        top.cell.set_box(frame.cell.vectors)
        top.set_positions(frame.positions)

    if args.ua:
        for mol in top.molecules:
            mol.remove_non_polar_hydrogens(update_topology=False)
        top.update_molecules(top.molecules)
        logger.info(str(top))

    if args.ff:
        ff = ForceField.open(args.ff)
        ff.assign_charge(top)
        logger.info(f'Net charge = {sum(a.charge for a in top.atoms)}')

    if args.qscale != 1:
        for mol in top.molecules:
            if mol.name in args.qscaleignore:
                continue
            for atom in mol.atoms:
                atom.charge *= args.qscale

    if args.box is not None:
        if len(args.box) == 3:
            box = args.box
        elif len(args.box) == 1:
            box = args.box * 3
        else:
            raise Exception('box should be float or list of three floats')
        top.cell.set_box(box)

    if top.has_position and set(args.shift) != {0}:
        for atom in top.atoms:
            atom.position += args.shift

    if not top.has_position:
        top.set_positions(np.zeros((top.n_atom, 3), dtype=float))

    top.write(args.output)
