#!/usr/bin/env python3

import os
import platform
import sys
import argparse
import math
import numpy as np
import matplotlib

if platform.system() == 'Linux':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})
from mstools.topology import Topology
from mstools.trajectory import Trajectory

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('cmd', choices=['dist', 'charge'],
                    help='the property to analyze')
parser.add_argument('-t', '--topology', required=True, type=str,
                    help='psf or lammps data file for topology information')
parser.add_argument('-i', '--input', nargs='+', required=True, type=str,
                    help='trajectory files for atomic positions and charges')
parser.add_argument('-b', '--begin', default=0, type=int,
                    help='first frame to analyze. Index starts from 0')
parser.add_argument('-e', '--end', default=-1, type=int,
                    help='last frame (not included) to analyze. Index starts from 0. '
                         '-1 means until last frames (included)')
parser.add_argument('--types', nargs='+', default=[], type=str,
                    help='consider these atom types')
parser.add_argument('--skip', default=1, type=int, help='skip frames in trajectory')
parser.add_argument('--dr', default=0.02, type=float, help='interval for density distribution')

args = parser.parse_args()

top = Topology.open(args.topology)
charges = np.array([atom.charge for atom in top.atoms], dtype=np.float32)

idx_zn1 = np.array([atom.id for atom in top.atoms if atom.symbol == 'Zn' and atom.molecule.id == 0])
idx_zn2 = np.array([atom.id for atom in top.atoms if atom.symbol == 'Zn' and atom.molecule.id == 1])
idx_p = np.array([atom.id for atom in top.atoms if atom.symbol == 'P'])
idx_n = np.array([atom.id for atom in top.atoms if atom.symbol == 'N' and atom.molecule.name == 'NTF2'])
idx_o = np.array([atom.id for atom in top.atoms if atom.symbol == 'O' and atom.molecule.name == 'PEG500'])
idx_list = [idx_p, idx_n, idx_o]

print(len(idx_zn1), len(idx_zn2))
print([len(l) for l in idx_list])

trj = Trajectory.open(args.input)

if args.end > trj.n_frame or args.end == -1:
    args.end = trj.n_frame


def calc_rsq(xyz1, xyz2):
    return np.sum((xyz2 - xyz1) ** 2)


def calc_r_z(vec, vec_unit):
    z = vec.dot(vec_unit)
    rsq = vec.dot(vec) - z ** 2
    if rsq < 0 and rsq > -1E6:
        r = 0
    else:
        r = np.sqrt(rsq)
    return r, z


dr = args.dr
max_r = 3
max_z = 5
min_z = -5
r_list = [(i + 0.5) * dr for i in range(int(max_r / dr))]
z_list = [(i + 0.5) * dr + min_z for i in range(int((max_z - min_z) / dr))]
densities = [np.zeros([len(z_list), len(r_list)]) for l in idx_list]
density_charge = np.zeros([len(z_list), len(r_list)])

points = []

n_frame = 0
for i in range(args.begin, args.end, args.skip):
    sys.stdout.write('\r\t %i / %i' % (i + 1, trj.n_frame))

    frame = trj.read_frame(i)
    n_frame += 1

    pos_zn1 = frame.positions[idx_zn1]
    com_zn1 = np.mean(pos_zn1, axis=0)

    pos_zn2 = frame.positions[idx_zn2]
    com_zn2 = np.mean(pos_zn2, axis=0)

    com_zn = np.mean([com_zn1, com_zn2], axis=0)
    vec_com = com_zn2 - com_zn
    vec_unit = vec_com / np.sqrt(vec_com.dot(vec_com))

    # print(com_zn1, com_zn2, vec, vec_unit)

    # consider PBC
    for mol in top.molecules:
        pos_mol = frame.positions[[atom.id for atom in mol.atoms]]
        com_mol = np.mean(pos_mol, axis=0)
        displace = com_mol - com_zn
        delta = np.floor(displace / frame.cell.lengths + 0.5)
        frame.positions[[atom.id for atom in mol.atoms]] -= delta * frame.cell.lengths

    if args.cmd == 'dist':
        for idx, density in zip(idx_list, densities):
            rho = len(idx) / frame.cell.volume
            for i in idx:
                vec = frame.positions[i] - com_zn
                _r, _z = calc_r_z(vec, vec_unit)
                if _r < max_r and _z < max_z and _z > min_z:
                    idx_r = int(_r / dr)
                    density[int((_z - min_z) / dr)][idx_r] += 1 / (2 * np.pi * r_list[idx_r] * dr * dr) / rho

    elif args.cmd == 'charge':
        for i in range(top.n_atom):
            vec = frame.positions[i] - com_zn
            _r, _z = calc_r_z(vec, vec_unit)
            if _r < max_r and _z < max_z and _z > min_z:
                idx_r = int(_r / dr)
                density_charge[int((_z - min_z) / dr)][idx_r] += charges[i] / (2 * np.pi * r_list[idx_r] * dr * dr)

if args.cmd == 'dist':
    for i, density in enumerate(densities):
        fig, ax = plt.subplots(1, 1, figsize=(5, 10))
        ax.set(xlim=[0, max_r], ylim=[min_z, max_z])
        im = ax.pcolormesh(r_list, z_list, density / n_frame, shading='nearest', vmin=0, vmax=4)
        fig.colorbar(im, ax=ax)
        fig.savefig('dist-%i.png' % i)
        if platform.system() != 'Linux':
            plt.show()

elif args.cmd == 'charge':
    fig, ax = plt.subplots(1, 1, figsize=(5, 10))
    ax.set(xlim=[0, max_r], ylim=[min_z, max_z])
    im = ax.pcolormesh(r_list, z_list, density_charge / n_frame, shading='nearest', cmap='PiYG', vmin=-20, vmax=20)
    fig.colorbar(im, ax=ax)
    fig.savefig('dist-charge.png')
    if platform.system() != 'Linux':
        plt.show()
