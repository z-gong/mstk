#!/usr/bin/env python3

import sys
import argparse
import math
import configparser
import numpy as np
import scipy.integrate as itg

np.seterr(all='raise')
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': 'Serif', 'font.size': 15})

from mstools.utils import histogram, print_data_to_file
from mstools.topology import Atom, Molecule, Topology
from mstools.trajectory import Trajectory
from mstools.constant import VACUUM_PERMITTIVITY, ELEMENTARY_CHARGE, NANO

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('cmd', choices=['density', 'diffuse'],
                    help='the property to analyze')
parser.add_argument('-t', '--topology', required=True, type=str,
                    help='psf or lammps data file for topology information')
parser.add_argument('-i', '--input', nargs='+', required=True, type=str,
                    help='trajectory files for atomic positions and charges')
parser.add_argument('-c', '--config', type=str, help='Config file for analysis')
parser.add_argument('-o', '--output', required=True, type=str, help='Output prefix')
parser.add_argument('-b', '--begin', default=0, type=int,
                    help='first frame to analyze. Index starts from 0')
parser.add_argument('-e', '--end', default=-1, type=int,
                    help='last frame (not included) to analyze. Index starts from 0. '
                         '-1 means until last frames (included)')
parser.add_argument('--dz', default=0.01, type=float,
                    help='bin size (nm) for numerical calculation of distribution and voltage')
parser.add_argument('--topignore', nargs='+', default=[], type=str,
                    help='ignore these molecule types in topology in case topology and trajectory do not match')
parser.add_argument('--dt', default=10, type=float,
                    help='time interval (ps) between frames if not presented in trajectory. '
                         'Required for diffusion analysis')
parser.add_argument('--skip', default=1, type=int, help='skip frames in trajectory')
parser.add_argument('--atoms', nargs='+', default=[], type=str,
                    help='analyze these atom types')

args = parser.parse_args()

top = Topology.open(args.topology)
if args.topignore != []:
    molecules = [mol for mol in top.molecules if mol.name not in args.topignore]
    top.init_from_molecules(molecules)
print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')

if args.cmd in ():
    ini = configparser.ConfigParser()
    if args.config is None or ini.read(args.config) == []:
        raise Exception('Config file is required')
    mol_names = ini['molecules']['name'].split()
    if set(mol_names) - set((mol.name for mol in top.molecules)) != set():
        raise Exception('Not all molecules listed in config file found in topology')

trj = Trajectory.open(args.input)
print('Trajectory info: ', trj.n_atom, 'atoms;', trj.n_frame, 'frames')

if (top.n_atom != trj.n_atom):
    raise Exception('Number of atoms in topology and trajectory files do not match')

# The box is periodic, all bins must have the same interval so that density profile makes sense
frame0 = trj.read_frame(0)
area = frame0.cell.size[0] * frame0.cell.size[1]
lz = frame0.cell.size[2]
dz = args.dz
n_bin = int(round((lz / dz)))
dz = lz / n_bin
edges = np.arange(n_bin + 1, dtype=np.float32) * dz
z_array = (edges[1:] + edges[:-1]) / 2

if args.end > trj.n_frame or args.end == -1:
    args.end = trj.n_frame


def replace_type_label(type):
    d = {'NC': 'Im$^+$',
         'PH': 'I$^-$',
         'CT': 'CH$_3$'
         }
    label = d.get(type)
    return label or type


def density():
    zcorr_dict = {atype: [] for atype in args.atoms}
    atom_ids = {atype: [atom.id for atom in top.atoms if atom.type == atype]
                for atype in args.atoms}
    n_frame = 0
    for i in range(args.begin, args.end, args.skip):
        n_frame += 1
        frame = trj.read_frame(i)
        sys.stdout.write('\r    frame %i' % i)
        for atype, ids in atom_ids.items():
            zs = frame.positions[ids][:, 2]
            zs -= np.floor(zs / lz) * lz
            zcorr_dict[atype] += list(zs)

    fig, ax = plt.subplots()
    ax.set(xlim=(0, lz))
    for atype, l in zcorr_dict.items():
        x, y = histogram(l, edges)
        y = y / n_frame / area / dz
        ax.plot(x, y, '.', label=replace_type_label(atype))
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'dens-{args.output}.png')


if __name__ == '__main__':
    if args.cmd == 'density':
        density()
