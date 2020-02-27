#!/usr/bin/env python3

import sys
import argparse
import math
import numpy as np
import configparser

import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})

from mstools.utils import histogram
from mstools.trajectory.topology import Atom, Molecule
from mstools.trajectory import Topology, Trajectory

parser = argparse.ArgumentParser()
parser.add_argument('cmd', choices=['dist', 'diffuse', 'voltage', 'charge2d', 'charge3d'],
                    help='The property to analyze')
parser.add_argument('-t', '--topology', required=True, type=str,
                    help='psf or lammps data file for topology information')
parser.add_argument('-i', '--input', required=True, type=str,
                    help='gro or lammpstrj file for atomic positions and charges')
parser.add_argument('-c', '--config', required=True, type=str, help='Config file for analysis')
parser.add_argument('-o', '--output', required=True, type=str, help='Output prefix')
parser.add_argument('-b', '--begin', default=0, type=int, help='first frame to output')
parser.add_argument('-e', '--end', default=-1, type=int, help='last frame to output')
parser.add_argument('--skip', default=1, type=int, help='skip frames between output')
parser.add_argument('--ignore', default='', type=str, help='ignore these atom types')
parser.add_argument('--voltage', default=0, type=float,
                    help='voltage drop in 3d image charge simulation.' 'Required for charge3d analysis')
args = parser.parse_args()

eps0 = 8.854188E-12
q0 = 1.602176E-19
Ang = 1E-10

top = Topology.open(args.topology)
print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')
trj = Trajectory.open(args.input)
print('Trajectory info: ', trj.n_atom, 'atoms;', trj.n_frame, 'frames')

if (top.n_atom != trj.n_atom):
    raise Exception('Number of atoms in topology and trajectory files do not match')

ignore_list = args.ignore.split(',')

if args.end > trj.n_frame or args.end == -1:
    args.end = trj.n_frame

config = configparser.ConfigParser()
config.read(args.config)
molecules = config['molecules']['name'].split()


def _get_atoms(mol: Molecule, names: [str]):
    return [list(filter(lambda x: x.name == name, mol.atoms))[0] for name in names]


def _get_com_position(positions, atoms: [Atom]):
    atom_masses = np.array([atom.mass for atom in atoms])
    com_mass = sum(atom_masses)
    atom_ids = [atom.id for atom in atoms]
    atom_positions = positions[atom_ids]
    com_position = np.sum(atom_positions * np.transpose(np.array([atom_masses] * 3)), axis=0) / com_mass
    return com_position


def _calc_angle_xy(vec):
    '''
    Calculate the angle between a vector and xy plane
    '''
    vec_proj = np.array([vec[0], vec[1], 0])
    theta = np.arccos(vec.dot(vec_proj) / np.sqrt(vec.dot(vec)) / np.sqrt(vec_proj.dot(vec_proj))) * 180 / np.pi
    return theta


def _calc_acf(series):
    '''
    Calculate the auto correlation function of a series
    acf(t) = <h(0)h(t)>
    '''
    acf = []
    for delta in (range(int(len(series) * 0.5))):
        _tmp = []
        for t0 in range(len(series) - delta):
            v0 = series[t0]
            vt = series[t0 + delta]
            _tmp.append(v0 * vt)
        acf.append(np.mean(_tmp))

    return acf


def distribution():
    bins = np.linspace(0, 70, 701)
    dz = bins[1] - bins[0]
    z_array = (bins[1:] + bins[:-1]) / 2
    charge_array = np.array([0.] * len(z_array))
    z_atom_dict = {}
    z_com_dict = {}
    theta_dict = {}

    frame = trj.read_frame(0)
    area = frame.box[0] * frame.box[1]

    n_frame = 0
    for i in range(args.begin, args.end, args.skip):
        n_frame += 1
        frame = trj.read_frame(i)
        sys.stdout.write('\r    step %i' % frame.step)
        positions = frame.positions
        for mol in top.molecules:
            if mol.name not in molecules:
                continue

            for atom in mol.atoms:
                z = frame.positions[atom.id][2]
                q = frame.charges[atom.id] if frame.has_charge else atom.charge
                z_idx = math.floor((z - bins[0]) / dz)
                charge_array[z_idx] += q

            section = config['molecule.%s' % (mol.name)]
            dists = section['distributions'].split(';')
            for dist in dists:
                name, atoms = [x.strip() for x in dist.split(':')]
                if name not in z_atom_dict:
                    z_atom_dict[name] = []
                atoms = _get_atoms(mol, atoms.split())
                z_atom_dict[name] += [positions[atom.id][2] for atom in atoms]

            com_dists = section['com_distributions'].split(';')
            for dist in com_dists:
                name, atoms = [x.strip() for x in dist.split(':')]
                if name not in z_com_dict:
                    z_com_dict[name] = []
                atoms = _get_atoms(mol, atoms.split())
                com_pos = _get_com_position(positions, atoms)
                z_com_dict[name].append(com_pos[2])

            angles = section['angles'].split(';')
            for angle in angles:
                name, theta_atoms = [x.strip() for x in angle.split(':')]
                if name not in theta_dict:
                    theta_dict[name] = []
                com_atoms, z_range = [x.strip() for x in section['angle.%s.com_zrange' % name].split(':')]
                com_atoms = _get_atoms(mol, com_atoms.split())
                z_range = list(map(float, z_range.split()))
                com_pos = _get_com_position(positions, com_atoms)
                if com_pos[2] < z_range[0] or com_pos[2] > z_range[1]:
                    continue

                if len(theta_atoms.split()) == 2:
                    a1, a2 = _get_atoms(mol, theta_atoms.split())
                    vec = positions[a2.id] - positions[a1.id]
                    theta_dict[name].append(_calc_angle_xy(vec))
                elif len(theta_atoms.split()) == 3:
                    a1, a2, a3 = _get_atoms(mol, theta_atoms.split())
                    vec_12 = positions[a2.id] - positions[a1.id]
                    vec_13 = positions[a3.id] - positions[a1.id]
                    vec_normal = np.cross(vec_12, vec_13)
                    theta_dict[name].append(90 - _calc_angle_xy(vec_normal))
                else:
                    raise Exception('Invalid angle definition', name, theta_atoms)

    print('')

    fig, ax = plt.subplots()
    ax.set(xlim=[0, 70], xlabel='z (A)', ylabel='particle density (/$A^3$)')
    for name, z_list in z_atom_dict.items():
        x, y = histogram(z_list, bins=bins)
        ax.plot(x, y / area / dz / n_frame, label=name)
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-dist.png')

    fig, ax = plt.subplots()
    ax.set(xlim=[0, 70], xlabel='z (A)', ylabel='molecule density (/$A^3$)')
    ax2 = ax.twinx()
    ax2.set_ylabel('cumulative molecule number')
    for name, z_list in z_com_dict.items():
        x, y_com = histogram(z_list, bins=bins)
        ax.plot(x, y_com / area / dz / n_frame, label=name)
        ax2.plot(x, np.cumsum(y_com) / n_frame, '--', label=name)
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-dist-com.png')

    fig, ax = plt.subplots()
    ax.set(xlim=[0, 70], xlabel='z (A)', ylabel='charge density (e/$A^3$)')
    ax.plot(z_array, charge_array / area / dz / n_frame, label=' '.join(molecules))
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-charge.png')

    fig, ax = plt.subplots()
    ax.set(xlim=[0, 90], ylim=[0, 0.1], xlabel='theta', ylabel='probability')
    for name, t_list in theta_dict.items():
        x, y = histogram(t_list, bins=np.linspace(0, 90, 91), normed=True)
        ax.plot(x, y, label=name)
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-angle.png')


def diffusion():
    name_atoms_dict = {}  # {'ring': [[atom1, atom2, ...], [atom11, atom12, ...], ...]}
    t_list = []
    z_dict = {}  # {'ring': [[], [], ...]}
    residence_zrange_dict = {}
    residence_dict = {}  # {'ring': [[], [], ...]}
    acf_dict = {}  # {'ring': []}

    for mol in top.molecules:
        if mol.name not in molecules:
            continue

        diffusions = config['molecule.%s' % (mol.name)]['diffusions'].split(';')
        for diffusion in diffusions:
            name, atoms = [x.strip() for x in diffusion.split(':')]
            if name not in name_atoms_dict:
                name_atoms_dict[name] = []
                z_dict[name] = []
                residence_dict[name] = []
                residence_zrange_dict[name] = [float(x) for x in config['molecule.%s' % (mol.name)][
                    'diffusion.%s.residence_zrange' % name].split()]
            name_atoms_dict[name].append(_get_atoms(mol, atoms.split()))
            z_dict[name].append([])
            residence_dict[name].append([])

    for i in range(args.begin, args.end, args.skip):
        frame = trj.read_frame(i)
        sys.stdout.write('\r    step %i' % frame.step)

        t_list.append(frame.step / 1e6)  # ns
        for name, atoms_list in name_atoms_dict.items():
            for k, atoms in enumerate(atoms_list):
                com_position = _get_com_position(frame.positions, atoms)
                z_dict[name][k].append(com_position[2])
                residence = int(com_position[2] >= residence_zrange_dict[name][0]
                                and com_position[2] <= residence_zrange_dict[name][1])
                residence_dict[name][k].append(residence)

    for name, residence_series_list in residence_dict.items():
        acf_series_list = []
        for residence_series in residence_series_list:
            acf_series_list.append(_calc_acf(residence_series))
        acf_dict[name] = np.mean(acf_series_list, axis=0) / np.mean(residence_series_list)

    print('')

    t_array = np.array(t_list) - t_list[0]
    for name, z_series_list in z_dict.items():
        fig, ax = plt.subplots(figsize=(6.4, 12.8))
        ax.set(ylim=[0, 69.496], xlabel='time (ns)', ylabel='z (A)')
        for z_series in z_series_list:
            ax.plot(t_array, z_series)
        fig.tight_layout()
        fig.savefig(f'{args.output}-diffusion-{name}.png')

    fig, ax = plt.subplots()
    ax.set(ylim=[0, 1.2], xlabel='time (ns)', ylabel='residence auto correlation')
    for name, z_series_list in z_dict.items():
        _len = len(acf_dict[name])
        ax.plot(t_array[:_len], acf_dict[name], label=name)
    ax.plot([0, t_array[-1] * 0.75], [np.exp(-1), np.exp(-1)], '--', label='$e^{-1}$')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-residence.png')


def voltage():
    frame = trj.read_frame(0)
    area = frame.box[0] * frame.box[1]
    dz = 0.1
    n_bin = math.ceil(frame.box[2] / dz)
    bins = np.array([frame.zlo + dz * (i + 0.5) for i in range(n_bin)])
    charges = np.array([0.] * n_bin)
    voltage = np.array([0.] * n_bin)

    n_frame = 0
    for i in range(args.begin, args.end, args.skip):
        n_frame += 1
        frame = trj.read_frame(i)
        sys.stdout.write('\r    step %i' % frame.step)

        for k, atom in enumerate(top.atoms):
            if atom.type in ignore_list:
                continue

            z = frame.positions[k][2]
            i_bin = math.floor((z - frame.zlo) / dz)
            i_bin = min(i_bin, n_bin - 1)
            q = frame.charges[k] if frame.has_charge else atom.charge
            charges[i_bin] += q

    charges_cumulative = np.cumsum(charges) / n_frame
    charges /= area * dz * n_frame  # e/A^3

    for i in range(1, n_bin):
        s = 0
        for j in range(0, i + 1):
            s += dz * (dz * (i - j)) * charges[j]
        voltage[i] = -s / eps0 * q0 / Ang

    fig = plt.figure(figsize=[6.4, 12.8])
    ax1 = fig.add_subplot('311')
    ax1.set(xlabel='z (A)', ylabel='charge density (e/A$^3$)')
    ax1.plot(bins, charges)
    ax1.plot(bins, [0] * n_bin, '--')

    ax2 = fig.add_subplot('312')
    ax2.set(ylim=[-5, 5], xlabel='z (A)', ylabel='cumulative charges (e)')
    ax2.plot(bins, charges_cumulative)
    ax2.plot(bins, [0] * n_bin, '--')

    ax3 = fig.add_subplot('313')
    ax3.set(xlabel='z (A)', ylabel='voltage (V)')
    ax3.plot(bins, voltage)
    ax3.plot(bins, [0] * n_bin, '--')
    fig.tight_layout()
    fig.savefig(f'{args.output}-voltage.png')


def charge_2d():
    frame = trj.read_frame(0)
    if not frame.has_charge:
        raise Exception('charge_2d function requires charge information in trajectory')

    z_span = frame.box[2] / 3
    _conv = q0 / frame.box[0] / frame.box[1] / Ang ** 2 * 1000  # convert from charge (e) to charge density (mC/m^2)
    ids_cathode = [atom.id for atom in top.atoms if atom.molecule.id == 1]
    qtot_list = []
    for i in range(args.begin, args.end, args.skip):
        frame = trj.read_frame(i)
        qtot = sum(frame.charges[ids_cathode])
        qtot_list.append(qtot)
        print('%-6i %10.6f %10.6f' % (i, qtot * _conv, qtot / len(ids_cathode) * 3))

    print('\n%-6i %10.6f %10.6f %10.6f' % (
        len(qtot_list), np.mean(qtot_list) * _conv, np.mean(qtot_list) / len(ids_cathode) * 3,
        np.mean(qtot_list) * _conv / 1000 * z_span * Ang / eps0))


def charge_3d():
    frame = trj.read_frame(0)
    area = frame.box[0] * frame.box[1]
    z_span = frame.box[2] / 2
    _conv = q0 / area / Ang ** 2 * 1000  # convert from charge (e) to charge density (mC/m^2)
    ids_cathode = [atom.id for atom in top.atoms if atom.molecule.id == 1]
    qtot_list = []
    for i in range(args.begin, args.end, args.skip):
        frame = trj.read_frame(i)
        qtot = args.voltage * area / z_span * eps0 / q0 * Ang
        for ii, atom in enumerate(top.atoms):
            if atom.molecule.name == 'MoS2' or atom.type == 'IMG':
                continue
            z = frame.positions[ii][2]
            q = frame.charges[ii] if frame.has_charge else atom.charge
            qtot += q * z / z_span
        qtot_list.append(qtot)
        print('%-6i %10.6f %10.6f' % (i, qtot * _conv, qtot / len(ids_cathode) * 3))

    print('\n%-6i %10.6f %10.6f %10.6f' % (
        len(qtot_list), np.mean(qtot_list) * _conv, np.mean(qtot_list) * _conv / len(ids_cathode) * 3,
        np.mean(qtot_list) * _conv / 1000 * z_span * Ang / eps0))


if __name__ == '__main__':
    if args.cmd == 'dist':
        distribution()
    elif args.cmd == 'diffuse':
        diffusion()
    elif args.cmd == 'voltage':
        voltage()
    elif args.cmd == 'charge2d':
        charge_2d()
    elif args.cmd == 'charge3d':
        charge_3d()
