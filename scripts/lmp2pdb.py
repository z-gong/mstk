#!/usr/bin/env python3

import sys
import argparse

from mstools.trajectory.lammps import LammpsData, LammpsTrj
from mstools.trajectory.pdb import PDB

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', required=True, type=str, help='Lammps data file for topology information')
parser.add_argument('-i', '--input', required=True, type=str, help='Lammps dump trajectory file for atomic positions')
parser.add_argument('-o', '--output', required=True, type=str, help='Output PDB file')
parser.add_argument('-b', '--begin', default=0, type=int, help='first frame to output')
parser.add_argument('-e', '--end', default=-1, type=int, help='last frame to output')
parser.add_argument('--skip', default=1, type=int, help='skip frames between output')
args = parser.parse_args()

top = LammpsData(args.DATA)
trj = LammpsTrj(args.INPUT)
pdb = PDB(args.OUTPUT, 'w')

print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')
print('Trajectory info: ', trj.n_atom, 'atoms;', trj.n_frame, 'frames')

if (top.n_atom != trj.n_atom):
    raise Exception('Number of atoms in topology and trajectory files do not match')

if args.end > trj.n_frame:
    args.end = trj.n_frame

for i in range(args.BEGIN, args.END, args.SKIP):
    sys.stdout.write('\r    %i' % i)
    frame = trj.read_frame(i)
    pdb.write_frame(top, frame)
pdb.close()
