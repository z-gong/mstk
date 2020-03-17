#!/usr/bin/env python3

import sys
import argparse
import numpy as np

from mstools.topology import Topology
from mstools.trajectory import Trajectory

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--topology', required=True, type=str,
                    help='psf or lammps data file for topology information')
parser.add_argument('-i', '--input', required=True, type=str, help='gro or lammpstrj file for atomic positions')
parser.add_argument('-o', '--output', required=True, type=str, help='output gro or pdb file')
parser.add_argument('-b', '--begin', default=0, type=int, help='first frame to output')
parser.add_argument('-e', '--end', default=-1, type=int, help='last frame to output')
parser.add_argument('--skip', default=1, type=int, help='skip frames between output')
parser.add_argument('--topignore', default='', type=str,
                    help='ignore these molecule types in topology in case topology and trajectory do not match')
parser.add_argument('--ignore', default='', type=str, help='ignore these molecule types')
parser.add_argument('--box', nargs=3, default=[-1, -1, -1], type=float, help='overwrite the box dimensions')
parser.add_argument('--shift', nargs=3, default=[0, 0, 0], type=float, help='shift the positions of all atoms')
args = parser.parse_args()

_top = Topology.open(args.topology)
top_ignore_list = args.topignore.split(',')
top = Topology()
top.init_from_molecules([mol for mol in _top.molecules if mol.name not in top_ignore_list])
print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')
trj = Trajectory.open(args.input)
print('Trajectory info: ', trj.n_atom, 'atoms;', trj.n_frame, 'frames')

if (top.n_atom != trj.n_atom):
    raise Exception('Number of atoms in topology and trajectory files do not match')

trj_out = Trajectory.open(args.output, 'w')

ignore_list = args.ignore.split(',')
if ignore_list != []:
    subset = [atom.id for atom in top.atoms if atom.molecule.name not in ignore_list]
else:
    subset = None

if args.end > trj.n_frame or args.end == -1:
    args.end = trj.n_frame

for i in range(args.begin, args.end, args.skip):
    sys.stdout.write('\r    %i' % i)
    frame = trj.read_frame(i)
    box = np.array([args.box[k] if args.box[k] != -1 else frame.box[k] for k in range(3)])
    frame.box = box
    if args.shift[0] != 0 or args.shift[1] != 0 or args.shift[2] != 0:
        frame.positions += np.array([args.shift] * top.n_atom)
    trj_out.write_frame(top, frame, subset=subset)
trj_out.close()
