#!/usr/bin/env python3

import sys
import argparse
import numpy as np

from mstools.topology import Topology
from mstools.trajectory import Trajectory

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='topology file')
parser.add_argument('-o', '--output', required=True, type=str, help='output trajectory file')
parser.add_argument('--ignore', nargs='+', default=[], type=str, help='ignore these molecule types')
parser.add_argument('--box', nargs=3, default=[-1, -1, -1], type=float,
                    help='overwrite the box dimensions')
parser.add_argument('--shift', nargs=3, default=[0, 0, 0], type=float,
                    help='shift the positions of all atoms')
args = parser.parse_args()

top = Topology.open(args.input)
if args.ignore != []:
    molecules = [mol for mol in top.molecules if mol.name not in args.ignore]
    top = Topology(molecules)
print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')

if args.shift[0] != 0 or args.shift[1] != 0 or args.shift[2] != 0:
    pos_shift = np.array([args.shift] * top.n_atom)
else:
    pos_shift = None

box = np.array([args.box[k] if args.box[k] != -1 else top.cell.size[k] for k in range(3)])
top.cell.set_box(box)
if pos_shift is not None:
    for atom in top.atoms:
        atom.position += pos_shift
top.write(args.output)
