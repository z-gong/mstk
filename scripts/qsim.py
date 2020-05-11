#!/usr/bin/env python3

import sys
import argparse
import numpy as np

from mstools.topology import Topology
from mstools.trajectory import Trajectory
from mstools.forcefield import FFSet
from mstools.simsys import System

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', nargs='+', required=True, type=str,
                    help='topology files for molecules (should contain positions)')
parser.add_argument('-n', '--number', nargs='+', type=int, help='number of molecules')
parser.add_argument('-t', '--trajectory', type=str, help='trajectory file for positions')
parser.add_argument('-f', '--forcefield', nargs='+', required=True, type=str,
                    help='forcefield files')
parser.add_argument('-o', '--output', required=True, type=str,
                    choices=['packmol', 'openmm', 'gromacs', 'lammps'],
                    help='output format')
parser.add_argument('--box', nargs=3, required=True, type=float,
                    help='periodic box size')
args = parser.parse_args()

if args.number == []:
    args.number = [1] * len(args.input)

molecules = []
for inp, n in zip(args.input, args.number):
    top = Topology.open(inp)
    molecules += top.molecules * n
top = Topology(molecules)
top.cell.set_box(args.box)

ff = FFSet.open(args.forcefield)

top.assign_mass_from_ff(ff)
top.assign_charge_from_ff(ff)

trj = Trajectory.open(args.trajectory)
frame = trj.read_frame(trj.n_frame - 1)
top.set_positions(frame.positions)
if frame.cell.volume != 0:
    top.cell.set_box(frame.cell.vectors)

system = System(top, ff)
if args.output == 'packmol':
    mol_numbers = top.get_unique_molecules()
    top.update_molecules(list(mol_numbers.keys()))
    top.scale_with_packmol(list(mol_numbers.values()))
elif args.output == 'gromacs':
    system.export_gmx()
