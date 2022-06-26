#!/usr/bin/env python3

import sys
import argparse
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mstk.topology import Topology
from mstk.trajectory import Trajectory
from mstk.utils import print_data_to_file
from mstk.topology.geometry import periodic_distance

np.seterr(all='raise')
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 15})

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', nargs='+', type=str,
                    help='trajectory files for atomic positions and charges')
parser.add_argument('-t', '--topology', required=True, type=str,
                    help='psf or lammps data file for topology information')
parser.add_argument('--a1', type=str, nargs='+', help='type of central atoms')
parser.add_argument('--a2', type=str, nargs='+', help='type of coordinate atoms')
parser.add_argument('--m1', type=str, help='name of central molecules')
parser.add_argument('--m2', type=str, help='name of coordinate molecules')
parser.add_argument('-o', '--output', required=True, type=str, help='Output prefix')
parser.add_argument('-b', '--begin', default=0, type=int,
                    help='first frame to analyze. Index starts from 0')
parser.add_argument('-e', '--end', default=-1, type=int,
                    help='last frame (not included) to analyze. Index starts from 0. '
                         '-1 means until last frames (included)')
parser.add_argument('--maxr', default=1.0, type=float,
                    help='max distance (nm) for calculation of distribution')
parser.add_argument('--dr', default=0.01, type=float,
                    help='bin size (nm) for calculation of distribution')
parser.add_argument('--ignore', nargs='+', default=[], type=str,
                    help='ignore these molecules in topology in case topology and trajectory do not match')
parser.add_argument('--skip', default=1, type=int, help='skip frames in trajectory')
args = parser.parse_args()

top = Topology.open(args.topology)
if args.ignore != []:
    molecules = [mol for mol in top.molecules if mol.name not in args.ignore]
    top.update_molecules(molecules)
print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')

trj = Trajectory.open(args.input)
print('Trajectory info: ', trj.n_atom, 'atoms;', trj.n_frame, 'frames')

if (top.n_atom != trj.n_atom):
    raise Exception('Number of atoms in topology and trajectory files do not match')

if args.m1 is None:
    group1 = [[atom] for atom in top.atoms if atom.type in args.a1]
else:
    group1 = [[atom for atom in mol.atoms if atom.type in args.a1]
              for mol in top.molecules if mol.name == args.m1]
ids_group1 = [[atom.id for atom in atoms] for atoms in group1]
masses_group1 = [[atom.mass for atom in atoms] for atoms in group1]

if args.m2 is None:
    group2 = [[atom] for atom in top.atoms if atom.type in args.a2]
else:
    group2 = [[atom for atom in mol.atoms if atom.type in args.a2]
              for mol in top.molecules if mol.name == args.m2]
ids_group2 = [[atom.id for atom in atoms] for atoms in group2]
masses_group2 = [[atom.mass for atom in atoms] for atoms in group2]

if args.end > trj.n_frame or args.end == -1:
    args.end = trj.n_frame

dr = args.dr
n_bin = math.ceil(args.maxr / dr) + 1
edges = np.array([dr * (i - 0.5) for i in range(n_bin + 1)], dtype=np.float32)
r_array = (edges[1:] + edges[:-1]) / 2
rdf_array = np.zeros(n_bin, dtype=np.float32)


def _get_weighted_center(positions, weight):
    return np.sum(positions * np.array(weight)[:, np.newaxis], axis=0) / sum(weight)


def _get_com_group(positions, ids_group, masses_group):
    com_group = []
    for ids, masses in zip(ids_group, masses_group):
        poss = positions[ids]
        if len(ids) == 1:
            pos = poss[0]
        else:
            pos = _get_weighted_center(poss, masses)
        com_group.append(pos)
    return com_group


n_frame = 0
for i in range(args.begin, args.end, args.skip):
    n_frame += 1
    frame = trj.read_frame(i)
    sys.stdout.write('\r    frame %i' % i)

    com_group1 = _get_com_group(frame.positions, ids_group1, masses_group1)
    com_group2 = _get_com_group(frame.positions, ids_group2, masses_group2)

    vol = frame.cell.volume
    density_pair = len(group1) * len(group2) / vol

    density_array = np.zeros(len(r_array), dtype=np.float32)
    for com1 in com_group1:
        for com2 in com_group2:
            distance = periodic_distance(com1, com2, frame.cell.size, args.maxr + dr / 2)
            if distance is not None and distance < args.maxr + dr / 2:
                idx_r = int((distance - edges[0]) / dr)
                density_array[idx_r] += 1

    rdf_array[1:] += density_array[1:] / (4 * np.pi * r_array[1:] ** 2 * dr) / density_pair

rdf_array = rdf_array / n_frame

name_column_dict = {'r'  : r_array,
                    'rdf': rdf_array}
print_data_to_file(name_column_dict, f'{args.output}-rdf.txt')

fig, ax = plt.subplots()
ax.set(xlabel='r (nm)', ylabel='RDF')
ax.plot(r_array, rdf_array, label='RDF')
ax.legend()
fig.tight_layout()
fig.savefig(f'{args.output}-rdf.png')
