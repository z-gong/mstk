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
from mstk import logger

np.seterr(all='raise')
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 15})


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--top', type=str, required=True, help='topology file')
    parser.add_argument('-c', '--conf', nargs='+', type=str, required=True, help='trajectory files')
    parser.add_argument('-o', '--output', default='rdf', type=str, help='Basename for output files')
    parser.add_argument('--a1', type=str, nargs='+', help='type of central atoms')
    parser.add_argument('--a2', type=str, nargs='+', help='type of coordinate atoms')
    parser.add_argument('--m1', type=str, help='name of central molecules')
    parser.add_argument('--m2', type=str, help='name of coordinate molecules')
    parser.add_argument('--type', default='all', type=str, choices=['all', 'inter', 'intra'],
                        help='calculate RDF between all pairs, inter-molecular pairs, or intra-molecular pairs')
    parser.add_argument('-b', '--begin', default=0, type=int,
                        help='first frame to analyze. Index starts from 0')
    parser.add_argument('-e', '--end', default=-1, type=int,
                        help='last frame (not included) to analyze. Index starts from 0. '
                             '-1 means until last frames (included)')
    parser.add_argument('--skip', default=1, type=int, help='read every N frames')
    parser.add_argument('--maxr', default=1.0, type=float,
                        help='max distance (nm) for calculation of distribution')
    parser.add_argument('--dr', default=0.01, type=float,
                        help='bin size (nm) for calculation of distribution')
    parser.add_argument('--ignore', nargs='+', default=[], type=str,
                        help='ignore these molecules in topology in case topology and trajectory do not match')

    return parser.parse_args()


def main():
    args = parse_args()
    top = Topology.open(args.top)
    if args.ignore:
        molecules = [mol for mol in top.molecules if mol.name not in args.ignore]
        top.update_molecules(molecules)
    logger.info(top)

    trj = Trajectory.open(args.conf)
    logger.info(trj)

    if top.n_atom != trj.n_atom:
        raise Exception('Number of atoms in topology and trajectory files do not match')

    if args.m1:
        group1 = [atom for mol in top.molecules if mol.name == args.m1 for atom in mol.atoms if atom.type in args.a1]
    else:
        group1 = [atom for atom in top.atoms if atom.type in args.a1]

    if args.m2:
        group2 = [atom for mol in top.molecules if mol.name == args.m2 for atom in mol.atoms if atom.type in args.a2]
    else:
        group2 = [atom for atom in top.atoms if atom.type in args.a2]

    if args.end > trj.n_frame or args.end == -1:
        args.end = trj.n_frame

    dr = args.dr
    n_bin = math.ceil(args.maxr / dr) + 1
    edges = np.array([dr * (i - 0.5) for i in range(n_bin + 1)], dtype=float)
    r_array = (edges[1:] + edges[:-1]) / 2
    rdf_array = np.zeros(n_bin, dtype=float)

    n_frame = 0
    for i in range(args.begin, args.end, args.skip):
        n_frame += 1
        frame = trj.read_frame(i)
        sys.stdout.write('\r    frame %i' % i)

        vol = frame.cell.volume
        box = frame.cell.get_size()
        density_pair = len(group1) * len(group2) / vol
        density_array = np.zeros(len(r_array), dtype=float)
        for atom1 in group1:
            pos1 = frame.positions[atom1.id]
            for atom2 in group2:
                if args.type == 'inter' and atom1.molecule is atom2.molecule:
                    continue
                if args.type == 'intra' and atom1.molecule is not atom2.molecule:
                    continue
                pos2 = frame.positions[atom2.id]
                distance = periodic_distance(pos1, pos2, box, args.maxr + dr / 2)
                if distance is not None and distance < args.maxr + dr / 2:
                    idx_r = int((distance - edges[0]) / dr)
                    density_array[idx_r] += 1

        rdf_array[1:] += density_array[1:] / (4 * np.pi * r_array[1:] ** 2 * dr) / density_pair

    rdf_array = rdf_array / n_frame

    name_column_dict = {'r'  : r_array,
                        'rdf': rdf_array}
    print_data_to_file(name_column_dict, f'{args.output}.txt')

    fig, ax = plt.subplots()
    ax.set(xlabel='r (nm)', ylabel='RDF')
    ax.plot(r_array, rdf_array)
    # ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}.png')


if __name__ == '__main__':
    main()
