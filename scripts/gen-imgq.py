#!/usr/bin/env python3

'''
Generate image charges for xyz file used for constant voltage simulation.
The z position of cathode and/or anode should be set.
The positions of images will be determined by treating electrodes as mirrors
'''

import sys
import argparse
from mstools.topology import XyzTopology
from mstools.trajectory import Xyz

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('input', required=True, type=str, help='input xyz file')
parser.add_argument('-o', '--output', required=True, type=str, help='output xyz file')
parser.add_argument('--cathode', type=float, help='z coordinate (in nm) of cathode (left electrode)')
parser.add_argument('--anode', type=float, help='z coordinate (in nm) of anode (right electrode)')
parser.add_argument('--ignore', default='', type=str, help='ignore these atom types')
parser.add_argument('--drude', action='store_true', help='generate image for drude particles for heavy atoms')
args = parser.parse_args()

top = XyzTopology(args.input)
trj = Xyz(args.input)
frame = trj.read_frame(0)

print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')
print('Trajectory info: ', trj.n_atom, 'atoms;', trj.n_frame, 'frames')

ignore_list = args.ignore.split(',')
gen_atoms = [atom for atom in top.atoms if atom.type not in ignore_list]
gen_drude_atoms = [atom for atom in gen_atoms if not atom.symbol == 'H'] if args.drude else []
for atom in top.atoms:
    atom._gen_img_ = atom._gen_img_drude_ = False
for atom in gen_atoms:
    atom._gen_img_ = True
for atom in gen_drude_atoms:
    atom._gen_img_drude_ = True

n_image = (len(gen_atoms) + len(gen_drude_atoms)) * ((args.cathode is not None) + (args.anode is not None))
if n_image == 0:
    print('ERROR: at least one of cathode and anode is required')
    sys.exit()
print('Generate %i image atoms' % n_image)

with open(args.output, 'w') as f_out:
    f_out.write('%i\n' % (top.n_atom + n_image))
    f_out.write('simulation box with image charges\n')
    for ii, atom in enumerate(top.atoms):
        pos = frame.positions[ii] * 10  # convert from nm to A
        f_out.write('%-8s %10.5f %10.5f %10.5f\n' % (atom.type, pos[0], pos[1], pos[2]))
    if args.cathode is not None:
        for ii, atom in enumerate(top.atoms):
            pos = frame.positions[ii] * 10  # convert from nm to A
            if atom._gen_img_:
                f_out.write('%-8s %10.5f %10.5f %10.5f\n' % ('IMG', pos[0], pos[1], args.cathode * 10 - pos[2]))
            if atom._gen_img_drude_:
                f_out.write('%-8s %10.5f %10.5f %10.5f\n' % ('IMG', pos[0], pos[1], args.cathode * 10 - pos[2]))
    if args.anode is not None:
        for ii, atom in enumerate(top.atoms):
            pos = frame.positions[ii] * 10  # convert from nm to A
            if atom._gen_img_:
                f_out.write('%-8s %10.5f %10.5f %10.5f\n' % ('IMG', pos[0], pos[1], 2 * args.anode * 10 - pos[2]))
            if atom._gen_img_drude_:
                f_out.write('%-8s %10.5f %10.5f %10.5f\n' % ('IMG', pos[0], pos[1], 2 * args.anode * 10 - pos[2]))
