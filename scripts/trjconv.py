#!/usr/bin/env python3

import sys
import argparse
import numpy as np

from mstk.topology import Topology
from mstk.trajectory import Trajectory
from mstk import logger


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs='+', type=str, help='trajectory files')
    parser.add_argument('-p', '--top', required=True, type=str, help='topology file')
    parser.add_argument('-o', '--output', required=True, type=str, help='output trajectory file')
    parser.add_argument('-b', '--begin', default=0, type=int, help='first frame to read')
    parser.add_argument('-e', '--end', default=-1, type=int, help='last frame (not included unless set to -1) to read')
    parser.add_argument('--skip', default=1, type=int, help='read every N frames')
    parser.add_argument('--topignore', nargs='+', default=[], type=str,
                        help='ignore these molecule names in topology in case topology and trajectory do not match')
    parser.add_argument('--ignore', nargs='+', default=[], type=str, help='ignore these molecule names')
    parser.add_argument('--ignoreatom', nargs='+', default=[], type=str, help='ignore these atom types')
    parser.add_argument('--box', nargs=3, default=[-1, -1, -1], type=float,
                        help='overwrite the box dimensions')
    parser.add_argument('--shift', nargs=3, default=[0, 0, 0], type=float,
                        help='shift the positions of all atoms')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    top = Topology.open(args.top)
    if args.topignore != []:
        molecules = [mol for mol in top.molecules if mol.name not in args.topignore]
        top = Topology()
        top.update_molecules(molecules)
    logger.info(str(top))

    trj = Trajectory.open(args.input)
    logger.info(str(trj))

    if (top.n_atom != trj.n_atom):
        raise Exception('Number of atoms in topology and trajectory files do not match')

    trj_out = Trajectory.open(args.output, 'w')

    if args.ignore or args.ignoreatom:
        subset = [atom.id for atom in top.atoms
                  if atom.type not in args.ignoreatom and atom.molecule.name not in args.ignore]
    else:
        subset = None

    if args.end > trj.n_frame or args.end == -1:
        args.end = trj.n_frame

    if args.shift[0] != 0 or args.shift[1] != 0 or args.shift[2] != 0:
        pos_shift = np.array([args.shift] * top.n_atom)
    else:
        pos_shift = None

    for i in range(args.begin, args.end, args.skip):
        sys.stdout.write('\r    %i' % i)
        frame = trj.read_frame(i)
        box = np.array([args.box[k] if args.box[k] != -1 else frame.cell.size[k] for k in range(3)])
        frame.cell.set_box(box)
        if pos_shift is not None:
            frame.positions += pos_shift
        trj_out.write_frame(frame, top, subset=subset)
    trj_out.close()
