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
parser.add_argument('--anode', required=True, type=float, help='location (z coordinate) of anode')
parser.add_argument('--cathode', type=float,help='location (z coordinate) of cathode if use Petersen\'s method')
parser.add_argument('--ignore', default='', type=str, help='ignore these atom types')
args = parser.parse_args()

top = XYZTopology(args.input)
trj = XYZ(args.input)
trj_out = XYZ(args.output, 'w')

print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')
print('Trajectory info: ', trj.n_atom, 'atoms;', trj.n_frame, 'frames')

ignore_list = args.ignore.split(',')
id_atoms = [atom.id for atom in top.atoms if atom.type not in ignore_list]

n_images = len(id_atoms)
if args.cathode is not None:
    n_images = 2 * len(id_atoms)

frame = trj.read_frame(0)
frame_out = Frame(top.n_atom + n_images)

for ii in range(top.n_atom):
    frame_out.positions[ii] = frame.positions[ii]

for ii, id in enumerate(id_atoms):
    atom = Atom()
    atom.type = 'IMG'
    top.atoms.append(atom)
    top.n_atom += 1
    frame_out.positions[ii + len(frame.positions)] = frame.positions[id]
    frame_out.positions[ii + len(frame.positions)][2] = 2 * args.anode - frame.positions[id][2]

if args.cathode is not None:
    for ii, id in enumerate(id_atoms):
        atom = Atom()
        atom.type = 'IMG'
        top.atoms.append(atom)
        top.n_atom += 1
        frame_out.positions[ii + len(frame.positions) + len(id_atoms)] = frame.positions[id]
        frame_out.positions[ii + len(frame.positions) + len(id_atoms)][2] = 2 * args.cathode - frame.positions[id][2]

trj_out.write_frame(top, frame_out)
trj_out.close()
