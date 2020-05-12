#!/usr/bin/env python3

import argparse
from mstools.topology import Topology
from mstools.trajectory import Trajectory
from mstools.constant import ELEMENTARY_CHARGE, NANO, VACUUM_PERMITTIVITY

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', type=str, help='input topology file')
parser.add_argument('-o', '--output', required=True, type=str, help='output topology file')
parser.add_argument('-q', '--charge', required=True, type=float, help='charge density (in mC/m^2)')
parser.add_argument('--trj', type=str, help='trajectory file for positions and cell')
parser.add_argument('--atom', default='SMo', type=str, help='type of electrode atoms')
parser.add_argument('--cathode', default=0.0, type=float, help='z coordinate (in nm) of cathode')
parser.add_argument('--anode', default=6.5, type=float, help='z coordinate (in nm) of anode')
parser.add_argument('--tolerance', default=0.05, type=float,
                    help='position tolerance for identifying electrode atoms')
args = parser.parse_args()

if args.output == args.input:
    raise Exception('Output filename should be different from input')

top = Topology.open(args.input)
if not top.has_position or top.cell.volume == 0:
    if args.trj is None:
        raise Exception('Position or box not found in topology file, trajectory should be provided')

    frame = Trajectory.read_frame_from_file(args.trj, 0)
    top.set_positions(frame.positions)
    top.cell.set_box(frame.cell.vectors)

atoms_catho = [atom for atom in top.atoms if atom.type == args.atom
               and abs(atom.position[2] - args.cathode) <= args.tolerance]
atoms_anode = [atom for atom in top.atoms if atom.type == args.atom
               and abs(atom.position[2] - args.anode) <= args.tolerance]

area = top.cell.size[0] * top.cell.size[1]
q_total = args.charge / 1000 * area * NANO ** 2 / ELEMENTARY_CHARGE
voltage = args.charge / 1000 * (args.anode - args.cathode) * NANO / VACUUM_PERMITTIVITY
q_catho = q_total / len(atoms_catho)
q_anode = -q_total / len(atoms_anode)

print('%.6f charges spread on each electrode' % q_total)
print('%.3f V voltage drop generated' % voltage)
print('%i atoms with extra charge %10.6f on cathode' % (len(atoms_catho), q_catho))
print('%i atoms with extra charge %10.6f on anode' % (len(atoms_anode), q_anode))

for atom in atoms_catho:
    atom.charge += q_catho
for atom in atoms_anode:
    atom.charge += q_anode

top.write(args.output)
