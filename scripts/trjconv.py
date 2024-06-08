#!/usr/bin/env python3

import sys
import argparse
import numpy as np

from mstk.topology import Topology
from mstk.trajectory import Trajectory
from mstk import logger


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--top', type=str, required=True, help='topology file')
    parser.add_argument('-c', '--conf', nargs='+', type=str, required=True, help='trajectory files')
    parser.add_argument('-o', '--output', required=True, type=str, help='output trajectory file')
    parser.add_argument('-b', '--begin', default=0, type=int,
                        help='read from this frame (included). '
                             'A negative value has the same meaning as in Python list')
    parser.add_argument('-e', '--end', type=int,
                        help='read until this frame (not included). '
                             'A negative value has the same meaning as in Python list. '
                             'If not set, will read until the end of the trajectory')
    parser.add_argument('--skip', default=1, type=int, help='read every N frames')
    parser.add_argument('--ignore', nargs='+', default=[], type=str, help='ignore these molecule names')
    parser.add_argument('--ignoreatom', nargs='+', default=[], type=str, help='ignore these atom types')
    parser.add_argument('--wrapfirst', action='store_true',
                        help='shift the positions of all atoms so that all residues in the first frame are in the main cell')
    parser.add_argument('--box', nargs=3, default=[-1, -1, -1], type=float,
                        help='overwrite the box dimensions')
    parser.add_argument('--shift', nargs=3, default=[0, 0, 0], type=float,
                        help='shift the positions of all atoms')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    top = Topology.open(args.top)
    logger.info(top)

    trj = Trajectory.open(args.conf)
    logger.info(trj)

    if (top.n_atom != trj.n_atom):
        raise Exception('Number of atoms in topology and trajectory files do not match')

    trj_out = Trajectory.open(args.output, 'w')

    if args.ignore or args.ignoreatom:
        subset = [atom.id for atom in top.atoms
                  if atom.type not in args.ignoreatom and atom.molecule.name not in args.ignore]
    else:
        subset = None

    if args.begin < 0:
        args.begin += trj.n_frame
    if args.end is None or args.end > trj.n_frame:
        args.end = trj.n_frame
    elif args.end < 0:
        args.end += trj.n_frame

    cell_offset = np.zeros((top.n_atom, 3), dtype=int)
    if args.wrapfirst:
        frame = trj.read_frame(args.begin)
        box = frame.cell.get_size()
        for res in top.residues:
            ids = [atom.id for atom in res.atoms]
            positions = frame.positions[ids]
            center = np.sum(positions, axis=0) / len(positions)
            cell_offset[ids] = np.floor(center / box).astype(int)

    if any(val != 0 for val in args.shift):
        pos_shift = np.array([args.shift] * top.n_atom)
    else:
        pos_shift = None

    for i in range(args.begin, args.end, args.skip):
        sys.stdout.write('\r    %i' % i)
        frame = trj.read_frame(i)
        box = frame.cell.get_size()
        for k in range(3):
            if args.box[k] != -1:
                box[k] = args.box[k]
        frame.cell.set_box(box)
        if args.wrapfirst:
            frame.positions -= cell_offset * box
        if pos_shift is not None:
            frame.positions += pos_shift
        trj_out.write_frame(frame, top, subset=subset)
    trj_out.close()
