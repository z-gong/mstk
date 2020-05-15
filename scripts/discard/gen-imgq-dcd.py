#!/usr/bin/env python3

'''
Generate image charges for xyz file used for constant voltage simulation.
The z position of cathode and/or anode should be set.
The positions of images will be determined by treating electrodes as mirrors
'''

import sys
import argparse
from mstools.topology import Topology
from mstools.trajectory import Trajectory

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-t', '--topology', required=True, type=str, help='topology file')
parser.add_argument('-i', '--input', required=True, type=str, help='input dcd file')
parser.add_argument('-o', '--output', required=True, type=str, help='output dcd file')
parser.add_argument('--anode', required=True, type=float,
                    help='z coordinate (in nm) of anode (right electrode)')
parser.add_argument('--ignore', nargs='+', default=[], type=str, help='ignore these atom types')
args = parser.parse_args()

top = Topology.open(args.topology)
trj = Trajectory.open(args.input)
frame = trj.read_frame(0)

atoms_non_img = [atom for atom in top.atoms if atom.type != 'IMG']
ids_gen = [atom.id for atom in atoms_non_img if atom.type not in args.ignore]

print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')
print('Trajectory info: ', trj.n_atom, 'atoms;', trj.n_frame, 'frames')
print('Image info: ', len(atoms_non_img), 'non-image atoms;', len(ids_gen), 'images to gen')

assert len(atoms_non_img) == trj.n_atom
assert trj.n_atom + len(ids_gen) == top.n_atom

trj_out = Trajectory.open(args.output, 'w')

for i in range(trj.n_frame):
    if i % 100 == 0:
        sys.stdout.write('\r    frame %i' % i)
    frame = trj.read_frame(i)
    frame.resize(top.n_atom)
    img_positions = frame.positions[ids_gen]
    img_positions[:, 2] = 2 * args.anode - img_positions[:, 2]
    frame.positions[trj.n_atom:] = img_positions
    trj_out.write_frame(frame, top)
