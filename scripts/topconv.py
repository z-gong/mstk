#!/usr/bin/env python3

import sys
import argparse
import numpy as np

from mstk.topology import Topology
from mstk.trajectory import Trajectory

parser = argparse.ArgumentParser()
parser.add_argument('input', nargs='+', type=str, help='topology file')
parser.add_argument('-n', '--number', nargs='+', type=int, help='number of molecules')
parser.add_argument('-o', '--output', required=True, type=str, help='output trajectory file')
parser.add_argument('--ignore', nargs='+', default=[], type=str, help='ignore these molecule types')
parser.add_argument('--qscale', default=1, type=float, help='scale the charge of atoms')
parser.add_argument('--qscaleignore', nargs='+', default=[], type=str,
                    help='ignore these molecule names for charge scaling')
parser.add_argument('--qscaleignoreatom', nargs='+', default=[], type=str,
                    help='ignore these atom types for charge scaling')
parser.add_argument('--box', nargs=3, default=[-1, -1, -1], type=float,
                    help='overwrite the box dimensions')
parser.add_argument('--shift', nargs=3, default=[0, 0, 0], type=float,
                    help='shift the positions of all atoms')
args = parser.parse_args()

top_list = [Topology.open(inp) for inp in args.input]
if args.number is None:
    args.number = [1] * len(top_list)

molecules = []
for top, n in zip(top_list, args.number):
    molecules.extend(top.molecules * n)
if args.ignore != []:
    molecules = [mol for mol in molecules if mol.name not in args.ignore]

top = top_list[0]
top.update_molecules(molecules)
print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')

if args.qscale != 1:
    for atom in top.atoms:
        if atom.type in args.qscaleignoreatom or atom.molecule.name in args.qscaleignore:
            continue
        atom.charge *= args.qscale

box = [args.box[k] if args.box[k] != -1 else top.cell.size[k] for k in range(3)]
top.cell.set_box(box)

if top.has_position and set(args.shift) != {0}:
    for atom in top.atoms:
        atom.position += args.shift

top.write(args.output)
