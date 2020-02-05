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
args = parser.parse_args()

top = XYZTopology(args.input)
trj = XYZ(args.input)
trj_out = XYZ(args.output, 'w')

print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')
print('Trajectory info: ', trj.n_atom, 'atoms;', trj.n_frame, 'frames')

ignore_list = args.ignore.split(',')
id_atoms = [atom.id for atom in top.atoms if atom.type not in ignore_list]

n_atom = top.n_atom
n_image = len(id_atoms) if args.anode is None else len(id_atoms) * 2

frame = trj.read_frame(0)
frame_out = Frame(n_atom + n_image)

for ii in range(n_atom):
    frame_out.positions[ii] = frame.positions[ii]

for ii, id in enumerate(id_atoms):
    img = Atom()
    img.type = 'IMG'
    top.atoms.append(img)
    top.n_atom += 1
    frame_out.positions[ii + n_atom] = frame.positions[id]
    frame_out.positions[ii + n_atom][2] = 2 * args.cathode - frame.positions[id][2]

if args.anode is not None:
    for ii, id in enumerate(id_atoms):
        img = Atom()
        img.type = 'IMG'
        top.atoms.append(img)
        top.n_atom += 1
        frame_out.positions[ii + n_atom + len(id_atoms)] = frame.positions[id]
        frame_out.positions[ii + n_atom + len(id_atoms)][2] = 2 * args.anode - frame.positions[id][2]

trj_out.write_frame(top, frame_out)
trj_out.close()
