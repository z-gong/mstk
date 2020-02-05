#!/usr/bin/env python3

'''
Generate image charge from xyz file
Used for constant voltage simulation
'''

import argparse
from mstools.trajectory.topology import Atom, Molecule
from mstools.trajectory.xyz import XYZTopology, XYZ
from mstools.trajectory import Frame

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, type=str, help='input xyz file')
parser.add_argument('-o', '--output', required=True, type=str, help='output xyz file')
parser.add_argument('--cathode', required=True, type=float, help='z coordinate of cathode (left electrode)')
parser.add_argument('--anode', type=float, help='z coordinate of anode (right electrode)')
parser.add_argument('--ignore', default='', type=str, help='ignore these atom types')
parser.add_argument('--drude', type=bool, action='store_true', help='generate image for drude particles for heavy atoms')
args = parser.parse_args()

top = XYZTopology(args.input)
trj = XYZ(args.input)
frame = trj.read_frame(0)

print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')
print('Trajectory info: ', trj.n_atom, 'atoms;', trj.n_frame, 'frames')

ignore_list = args.ignore.split(',')
gen_atoms = [atom for atom in top.atoms if atom.type not in ignore_list]
hvy_atoms = [atom for atom in gen_atoms if not atom.type.starswith('H')] if args.drude else []

n_atom = top.n_atom
n_image = len(gen_atoms) + len(hvy_atoms)
if args.anode is not None:
    n_image *= 2

with open(args.output, 'w') as f_out:
    f_out.write('%i\n' % (n_atom + n_image))
    f_out.write('simbox with image charges\n')
    for ii, atom in enumerate(top.atoms):
        pos = frame.positions
        f_out.write('%-8s %10.5f %10.5f %10.5f\n' % (atom.type, pos[0], pos[1], pos[2]))
        if atom in gen_atoms:
            f_out.write('%-8s %10.5f %10.5f %10.5f\n' % (atom.type, pos[0], pos[1], args.cathode - pos[2]))
        if atom in hvy_atoms:
            f_out.write('%-8s %10.5f %10.5f %10.5f\n' % (atom.type, pos[0], pos[1], args.cathode - pos[2]))
    if args.anode is not None:
        for ii, atom in enumerate(top.atoms):
            pos = frame.positions
            if atom in gen_atoms:
                f_out.write('%-8s %10.5f %10.5f %10.5f\n' % (atom.type, pos[0], pos[1], 2 * args.anode - pos[2]))
            if atom in hvy_atoms:
                f_out.write('%-8s %10.5f %10.5f %10.5f\n' % (atom.type, pos[0], pos[1], 2 * args.anode - pos[2]))
