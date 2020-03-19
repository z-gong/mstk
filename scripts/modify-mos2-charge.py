#!/usr/bin/env python3

import argparse
from mstools.topology import Topology
from mstools.trajectory import Trajectory

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, type=str, help='input topology file')
parser.add_argument('-o', '--output', required=True, type=str, help='output topology file')
parser.add_argument('-q', '--charge', required=True, type=float, help='charge density (in mC/m^2)')
parser.add_argument('--trj', type=str, help='trajectory file used to identify electrode atoms')
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

    trj = Trajectory.open(args.trj)
    frame = trj.read_frame(0)
    top.set_positions(frame.positions)
    top.cell = frame.cell

top_out = Topology.open(args.output, 'w')
top_out.init_from_topology(top)
atoms_cat = [atom for atom in top_out.atoms if atom.type == args.atom
             and abs(atom.position[2] - args.cathode) <= args.tolerance]
atoms_ano = [atom for atom in top_out.atoms if atom.type == args.atom
             and abs(atom.position[2] - args.anode) <= args.tolerance]

e0 = 1.60217662E-19
nm = 1E-9
area = top.cell.size[0] * top.cell.size[1]
q_cat = args.charge / 1000 * area * nm ** 2 / len(atoms_cat) / e0
q_ano = args.charge / 1000 * area * nm ** 2 / len(atoms_ano) / e0

print('%i atoms with extra charge %.6f on cathode')
print('%i atoms with extra charge %.6f on anode')

for atom in atoms_cat:
    print(atom, atom.position, atom.charge)
    atom.charge += q_cat
for atom in atoms_ano:
    print(atom, atom.position, atom.charge)
    atom.charge -= q_ano

top_out.write()
top_out.close()
