#!/usr/bin/env python3

'''
Generate image charges for xyz file used for constant voltage simulation.
The z position of cathode and/or anode should be set.
The positions of images will be determined by treating electrodes as mirrors
'''

import sys
import argparse
from mstk.topology import Topology

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('input', type=str, help='input xyz file')
parser.add_argument('-o', '--output', required=True, type=str, help='output xyz file')
parser.add_argument('--cathode', type=float, help='z coordinate (in nm) of cathode (left electrode)')
parser.add_argument('--anode', type=float, help='z coordinate (in nm) of anode (right electrode)')
parser.add_argument('--ignore', nargs='+', default=[], type=str, help='ignore these atom types')
parser.add_argument('--drude', action='store_true', help='generate image for drude particles for heavy atoms')
parser.add_argument('--vsite', nargs='+', default=[], help='generate image for virtual sites for specified atoms')
args = parser.parse_args()

top = Topology.open(args.input)
positions = top.positions[:]

print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')

gen_atoms = [atom for atom in top.atoms if atom.type not in args.ignore]
gen_drude_atoms = [atom for atom in gen_atoms if atom.symbol != 'H'] if args.drude else []
gen_vsite_atoms = [atom for atom in gen_atoms if atom.type in args.vsite]
for atom in top.atoms:
    atom._gen_img_ = atom._gen_img_drude_ = atom._gen_img_vsite_ = False
for atom in gen_atoms:
    atom._gen_img_ = True
for atom in gen_drude_atoms:
    atom._gen_img_drude_ = True
for atom in gen_vsite_atoms:
    atom._gen_img_vsite_ = True

n_image = (len(gen_atoms) + len(gen_drude_atoms) + len(gen_vsite_atoms)) * ((args.cathode is not None) + (args.anode is not None))
if n_image == 0:
    print('ERROR: at least one of cathode and anode is required')
    sys.exit()
print('Generate %i image atoms' % n_image)

with open(args.output, 'w') as f_out:
    f_out.write('%i\n' % (top.n_atom + n_image))
    f_out.write('simulation box with image charges\n')
    for ii, atom in enumerate(top.atoms):
        pos = positions[ii] * 10  # convert from nm to A
        f_out.write('%-8s %10.5f %10.5f %10.5f\n' % (atom.type, pos[0], pos[1], pos[2]))
    if args.cathode is not None:
        for ii, atom in enumerate(top.atoms):
            pos = positions[ii] * 10  # convert from nm to A
            if atom._gen_img_:
                f_out.write('%-8s %10.5f %10.5f %10.5f\n' % ('IMG', pos[0], pos[1], args.cathode * 10 - pos[2]))
            if atom._gen_img_drude_:
                f_out.write('%-8s %10.5f %10.5f %10.5f\n' % ('IMG', pos[0], pos[1], args.cathode * 10 - pos[2]))
            if atom._gen_img_vsite_:
                f_out.write('%-8s %10.5f %10.5f %10.5f\n' % ('IMG', pos[0], pos[1], args.cathode * 10 - pos[2]))
    if args.anode is not None:
        for ii, atom in enumerate(top.atoms):
            pos = positions[ii] * 10  # convert from nm to A
            if atom._gen_img_:
                f_out.write('%-8s %10.5f %10.5f %10.5f\n' % ('IMG', pos[0], pos[1], 2 * args.anode * 10 - pos[2]))
            if atom._gen_img_drude_:
                f_out.write('%-8s %10.5f %10.5f %10.5f\n' % ('IMG', pos[0], pos[1], 2 * args.anode * 10 - pos[2]))
            if atom._gen_img_vsite_:
                f_out.write('%-8s %10.5f %10.5f %10.5f\n' % ('IMG', pos[0], pos[1], 2 * args.anode * 10 - pos[2]))
