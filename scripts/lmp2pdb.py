#!/usr/bin/env python3

import sys
import argparse

from mstools.trajectory.lammps import LammpsData, LammpsTrj
from mstools.trajectory.pdb import PDB

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', required=True, type=str, help='lammps data file for topology information')
parser.add_argument('-i', '--input', required=True, type=str, help='lammps dump trajectory file for atomic positions')
parser.add_argument('-o', '--output', required=True, type=str, help='output PDB file')
parser.add_argument('-b', '--begin', default=0, type=int, help='first frame to output')
parser.add_argument('-e', '--end', default=-1, type=int, help='last frame to output')
parser.add_argument('--skip', default=1, type=int, help='skip frames between output')
parser.add_argument('--ignore', default='', type=str, help='ignore these atom types')
parser.add_argument('--box', nargs=3, default=[-1, -1, -1], type=float, help='overwrite the box dimensions')
args = parser.parse_args()

top = LammpsData(args.data)
print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')
trj = LammpsTrj(args.input)
print('Trajectory info: ', trj.n_atom, 'atoms;', trj.n_frame, 'frames')
pdb = PDB(args.output, 'w')

ignore_list = args.ignore.split(',')
id_atoms = [atom.id for atom in top.atoms if atom.type not in ignore_list]

if (top.n_atom != trj.n_atom):
    raise Exception('Number of atoms in topology and trajectory files do not match')

if args.end > trj.n_frame or args.end == -1:
    args.end = trj.n_frame

for i in range(args.begin, args.end, args.skip):
    sys.stdout.write('\r    %i' % i)
    frame = trj.read_frame(i)
    box = [args.box[k] if args.box[k] != -1 else frame.box[k] for k in range(3)]
    frame.box = box
    pdb.write_frame(top, frame, subset=id_atoms)
pdb.close()