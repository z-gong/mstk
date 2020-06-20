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

plt.rcParams.update({'font.size': 15})

from mstools.utils import histogram, print_data_to_file
from mstools.topology import Atom, Molecule, Topology
from mstools.trajectory import Trajectory
from mstools.constant import *

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('cmd', choices=['dist', 'diffuse', 'voltage', 'drude', 'dipole', 'charge',
                                    'fluctq', 'permit'],
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
parser.add_argument('--cathode', default=0, type=float, help='z coordinate (in nm) of cathode')
parser.add_argument('--anode', default=6.5, type=float, help='z coordinate (in nm) of anode')
parser.add_argument('--dz', default=0.01, type=float,
                    help='bin size (nm) for numerical calculation of distribution and voltage')
parser.add_argument('--ignore', nargs='+', default=[], type=str,
                    help='ignore these molecules in topology in case topology and trajectory do not match')
parser.add_argument('--dt', default=10, type=float,
                    help='time interval (ps) between frames if not presented in trajectory. '
                         'Required for diffusion analysis')
parser.add_argument('--skip', default=1, type=int, help='skip frames in trajectory')
parser.add_argument('--voltage', default=0, type=float,
                    help='pre-assigned voltage drop in image charge simulation for analyzing surface charge. '
                         'Required for charge analysis')
parser.add_argument('--charge', default=0, type=float,
                    help='pre-assigned or evolved charge density on electrodes for analyzing voltage profile. '
                         'Required for voltage analysis')
parser.add_argument('--reverse', action='store_true',
                    help='reverse the direction of cathode and anode. '
                         'It makes the analysis easier if anode is to be emphasized. '
                         'It only affect how the results are presented')

args = parser.parse_args()

top = Topology.open(args.topology)
if args.ignore != []:
    molecules = [mol for mol in top.molecules if mol.name not in args.ignore]
    top.update_molecules(molecules)
print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')

if args.cmd in ('dist', 'diffuse', 'dipole'):
    ini = configparser.ConfigParser()
    if args.config is None:
        raise Exception('Config file is required for dist, diffuse and dipole analysis')
    if ini.read(args.config) == []:
        raise Exception('Invalid config file')
    mol_names = ini['molecules']['name'].split()
    if set(mol_names) - set((mol.name for mol in top.molecules)) != set():
        raise Exception('Not all molecules listed in config file found in topology')

trj = Trajectory.open(args.input)
print('Trajectory info: ', trj.n_atom, 'atoms;', trj.n_frame, 'frames')

if (top.n_atom != trj.n_atom):
    raise Exception('Number of atoms in topology and trajectory files do not match')

frame0 = trj.read_frame(0)
ids_cathode = [i for i in range(len(frame0.positions)) if
               abs(frame0.positions[i][2] - args.cathode) < 0.05]
ids_anode = [i for i in range(len(frame0.positions)) if
             abs(frame0.positions[i][2] - args.anode) < 0.05]
if len(ids_cathode) == 0 or len(ids_anode) == 0 or len(ids_cathode) != len(ids_anode):
    print('WARNING: Cannot identify consistent atoms for cathode and anode. '
          'Make sure positions of cathode and anode are correct')
area = frame0.cell.size[0] * frame0.cell.size[1]

elecd_l = min(args.cathode, args.anode)
elecd_r = max(args.cathode, args.anode)
if args.reverse:
    cathode, anode = args.anode, args.cathode
else:
    cathode, anode = args.cathode, args.anode
lz = elecd_r - elecd_l
dz = args.dz
# increase both the bottom and top by 0.4 nm to consider electrodes
# add one extra bin and move edges by 0.5*dz so that the centers of bins are prettier
n_bin = math.ceil((lz + 0.8) / dz) + 1
edges = np.array([dz * (i - 0.5) + elecd_l - 0.4 for i in range(n_bin + 1)], dtype=np.float32)
z_array = (edges[1:] + edges[:-1]) / 2

if args.end > trj.n_frame or args.end == -1:
    args.end = trj.n_frame


def _get_atoms(mol: Molecule, names: [str]):
    return [next(atom for atom in mol.atoms if atom.name == name) for name in names]


def _get_weighted_center(positions, weight):
    return np.sum(positions * np.transpose(np.array([weight] * 3)), axis=0) / sum(weight)


def _get_com_position(positions, atoms: [Atom]):
    atom_ids = [atom.id for atom in atoms]
    atom_positions = positions[atom_ids]
    atom_masses = [atom.mass for atom in atoms]
    return _get_weighted_center(atom_positions, atom_masses)


def _calc_angle_xy(vec):
    '''
    Calculate the angle between a vector and xy plane
    '''
    if vec[0] == 0 and vec[1] == 0:
        theta = 90
    else:
        vec_prj = np.array([vec[0], vec[1], 0])
        cos = vec.dot(vec_prj) / np.sqrt(vec.dot(vec)) / np.sqrt(vec_prj.dot(vec_prj))
        theta = np.arccos(np.clip(cos, -1, 1)) * 180 / np.pi

    return theta


def _calc_acf(series):
    '''
    Calculate the auto correlation function of a series
    acf(t) = <(h(0)*h(t)>
    '''
    n = len(series)
    r = np.correlate(series, series, mode='full')[-n:]
    # assert np.allclose(r, np.array([(x[:n - k] * x[-(n - k):]).sum() for k in range(n)]))
    acf = r / np.arange(n, 0, -1)

    return acf


def distribution():
    z_atom_dict = {}
    z_com_dict = {}
    theta_dict = {}

    n_frame = 0
    for i in range(args.begin, args.end, args.skip):
        n_frame += 1
        frame = trj.read_frame(i)
        sys.stdout.write('\r    frame %i' % i)

        positions = frame.positions
        if args.reverse:
            positions[:, 2] = elecd_l + elecd_r - positions[:, 2]

        for mol in top.molecules:
            if mol.name not in mol_names:
                continue

            section = ini['molecule.%s' % (mol.name)]
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
                z_com_dict[name].append(_get_com_position(positions, atoms)[2])

            angles = section['angles'].split(';')
            for angle in angles:
                name, theta_atoms = [x.strip() for x in angle.split(':')]
                if name not in theta_dict:
                    theta_dict[name] = []
                com_atoms, z_range = [x.strip() for x in
                                      section['angle.%s.com_zrange' % name].split(':')]
                com_atoms = _get_atoms(mol, com_atoms.split())
                z_range = [float(x) + elecd_l for x in z_range.split()]
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

    name_column_dict = {'z': z_array}
    fig, ax = plt.subplots()
    ax.set(xlim=[edges[0], edges[-1]], xlabel='z (nm)', ylabel='particle density (/$nm^3$)')
    for name, z_list in z_atom_dict.items():
        x, y = histogram(z_list, bins=edges)
        ax.plot(x, y / area / dz / n_frame, label=name)
        name_column_dict['rho-' + name] = y / area / dz / n_frame
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-dist.png')
    print_data_to_file(name_column_dict, f'{args.output}-dist.txt')

    name_column_dict = {'z': z_array}
    fig, ax = plt.subplots()
    ax.set(xlim=[edges[0], edges[-1]], xlabel='z (nm)', ylabel='molecule density (/$nm^3$)')
    ax2 = ax.twinx()
    ax2.set_ylabel('cumulative molecule number')
    for name, z_list in z_com_dict.items():
        x, y_com = histogram(z_list, bins=edges)
        ax.plot(x, y_com / area / dz / n_frame, label=name,
                color='darkred' if name == 'dca_com' else None)
        ax2.plot(x, np.cumsum(y_com) / n_frame, '--', label=name,
                 color='darkred' if name == 'dca_com' else None)
        name_column_dict['rho-' + name] = y_com / area / dz / n_frame
        name_column_dict['cum-' + name] = np.cumsum(y_com) / n_frame
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-dist-com.png')
    print_data_to_file(name_column_dict, f'{args.output}-dist-com.txt')

    name_column_dict = {}
    fig, ax = plt.subplots()
    ax.set(xlim=[0, 90], ylim=[0, 0.1], xlabel='theta', ylabel='probability')
    for name, t_list in theta_dict.items():
        x, y = histogram(t_list, bins=np.linspace(0, 90, 91), normed=True)
        ax.plot(x, y, label=name, color='darkred' if name == 'dca_e2e' else None)
        name_column_dict.update({'theta': x, 'prob-' + name: y})
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-angle.png')
    print_data_to_file(name_column_dict, f'{args.output}-angle.txt')


def diffusion():
    name_atoms_dict = {}  # {'ring': [[atom1, atom2, ...], [atom11, atom12, ...], ...]}
    t_list = []
    z_dict = {}  # {'ring': [[], [], ...]}
    residence_zrange_dict = {}
    residence_dict = {}  # {'ring': [[], [], ...]}
    acf_dict = {}  # {'ring': []}

    for mol in top.molecules:
        if mol.name not in mol_names:
            continue

        diffusions = ini['molecule.%s' % (mol.name)]['diffusions'].split(';')
        for diffusion in diffusions:
            name, atoms = [x.strip() for x in diffusion.split(':')]
            if name not in name_atoms_dict:
                name_atoms_dict[name] = []
                z_dict[name] = []
                residence_dict[name] = []
                residence_zrange_dict[name] = [
                    float(x) + elecd_l for x in
                    ini['molecule.%s' % (mol.name)]['diffusion.%s.residence_zrange' % name].split()
                ]
            name_atoms_dict[name].append(_get_atoms(mol, atoms.split()))
            z_dict[name].append([])
            residence_dict[name].append([])

    for i in range(args.begin, args.end, args.skip):
        frame = trj.read_frame(i)
        sys.stdout.write('\r    frame %i' % i)

        positions = frame.positions
        if args.reverse:
            positions[:, 2] = elecd_l + elecd_r - positions[:, 2]

        if frame.time != -1:
            t_list.append(frame.time / 1E3)  # convert ps to ns
        else:
            t_list.append(i * args.dt / 1E3)

        for name, atoms_list in name_atoms_dict.items():
            for k, atoms in enumerate(atoms_list):
                com_position = _get_com_position(positions, atoms)
                z_dict[name][k].append(com_position[2])
                residence = int(com_position[2] >= residence_zrange_dict[name][0]
                                and com_position[2] <= residence_zrange_dict[name][1])
                residence_dict[name][k].append(residence)

    for name, residence_series_list in residence_dict.items():
        acf_series_list = []
        # for series in residence_series_list:
        #     acf_series_list.append(_calc_acf(series - np.mean(residence_series_list)))
        # acf_dict[name] = np.mean(acf_series_list, axis=0) / np.var(residence_series_list)
        for series in residence_series_list:
            acf_series_list.append(_calc_acf(series))
        acf_dict[name] = np.mean(acf_series_list, axis=0) / np.mean(residence_series_list)

    print('')

    t_array = np.array(t_list) - t_list[0]
    name_column_dict = {'time': t_array}
    for name, z_series_list in z_dict.items():
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 8))
        ax1.set(ylim=[elecd_r - 1.0, elecd_r], ylabel='z (nm)')
        ax2.set(ylim=[elecd_l, elecd_l + 1.0], xlabel='time (ns)', ylabel='z (nm)')
        ax1.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax1.xaxis.tick_top()
        ax2.xaxis.tick_bottom()
        for z_series in z_series_list:
            ax1.plot(t_array, z_series, linewidth=1)
            ax2.plot(t_array, z_series, linewidth=1)
        plt.subplots_adjust(hspace=0.1)
        fig.tight_layout()
        fig.savefig(f'{args.output}-diffusion-{name}.png')

    fig, ax = plt.subplots()
    ax.set(xlim=[0, t_array[-1] * 0.75], ylim=[0, 1.2], xlabel='time (ns)',
           ylabel='residence auto correlation')
    for name, z_series_list in z_dict.items():
        _len = len(acf_dict[name])
        ax.plot(t_array[:_len], acf_dict[name], label=name,
                color='darkred' if name == 'dca_com' else None)
        name_column_dict['time'] = t_array[:_len]
        name_column_dict['acf-' + name] = acf_dict[name]
    ax.plot([0, t_array[-1] * 0.75], [np.exp(-1), np.exp(-1)], '--', label='$e^{-1}$')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-residence.png')

    print_data_to_file(name_column_dict, f'{args.output}-residence.txt')


def drude_dipole():
    '''
    Calculate the distribution of Drude induced dipoles
    '''
    top_atom_charges = np.array([atom.charge for atom in top.atoms], dtype=np.float32)
    drude_pairs = top.get_drude_pairs()
    parent_ids = [parent.id for parent, drude in drude_pairs]
    drude_ids = [drude.id for parent, drude in drude_pairs]
    drude_charges = top_atom_charges[drude_ids]
    z_dipoles = [[] for z in z_array]
    z_dipoles_scalar = [[] for z in z_array]
    z_distances = [[] for z in z_array]
    n_frame = 0
    for i in range(args.begin, args.end, args.skip):
        n_frame += 1
        frame = trj.read_frame(i)
        sys.stdout.write('\r    frame %i' % i)

        positions = frame.positions
        if args.reverse:
            positions[:, 2] = elecd_l + elecd_r - positions[:, 2]

        drude_positions = positions[drude_ids]
        parent_positions = positions[parent_ids]
        deltas = drude_positions - parent_positions
        dipoles = deltas * drude_charges[:, np.newaxis]
        distances = np.array([np.sqrt(np.dot(d, d)) for d in deltas])
        dipoles_scalar = -distances * drude_charges[: np.newaxis]
        z_indexes = ((parent_positions[:, 2] - edges[0]) / dz).astype(int)
        for idx, dipole, dipole_scalar, dist in zip(z_indexes, dipoles, dipoles_scalar, distances):
            z_dipoles[idx].append(dipole)
            z_dipoles_scalar[idx].append(dipole_scalar)
            z_distances[idx].append(dist)

    z_dipole = np.array(
        [np.mean(dipoles, axis=0) if len(dipoles) > n_frame else np.array([0., 0., 0.])
         for dipoles in z_dipoles])
    z_distance = np.array([np.mean(distances) if len(distances) > n_frame else 0
                           for distances in z_distances])
    z_dipole_scalar = np.array([np.mean(dipoles_scalar) if len(dipoles_scalar) > n_frame else 0
                                for dipoles_scalar in z_dipoles_scalar])

    fig, ax = plt.subplots()
    ax.set(xlabel='z (nm)', ylabel='induced dipole (e*nm)')
    ax.plot(z_array, z_dipole_scalar, label='|dipole|')
    ax.plot(z_array, z_dipole[:, 2], label='dipoleZ')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-induced.png')

    fig, ax = plt.subplots()
    ax.set(xlabel='z (nm)', ylabel='Drude-parent distance (nm)')
    ax.plot(z_array, z_distance, label='distance')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-induced_distance.png')


def dipole():
    '''
    Calculate the distribution of molecular dipoles
    '''
    top_atom_charges = np.array([atom.charge for atom in top.atoms], dtype=np.float32)
    name_z_dipoles_dict = {}  # {'ring': [[array, array, ...], [array, array, ...] , ...]}
    name_z_dipole_dict = {}  # {'ring': array([array, array , ...])}
    n_frame = 0
    for i in range(args.begin, args.end, args.skip):
        n_frame += 1
        frame = trj.read_frame(i)
        sys.stdout.write('\r    frame %i' % i)

        positions = frame.positions
        if args.reverse:
            positions[:, 2] = elecd_l + elecd_r - positions[:, 2]

        for mol in top.molecules:
            if mol.name not in mol_names:
                continue

            section = ini['molecule.%s' % (mol.name)]
            dipoles = section['dipoles'].split(';')
            for d in dipoles:
                name, atoms = [x.strip() for x in d.split(':')]
                if name not in name_z_dipoles_dict:
                    name_z_dipoles_dict[name] = [[] for z in z_array]
                atoms = _get_atoms(mol, atoms.split())
                ids = [atom.id for atom in atoms]
                charges = frame.charges[ids] if frame.has_charge else top_atom_charges[ids]
                rel_charges = charges - charges.mean()
                poss = positions[ids]
                dipole = np.sum(poss * np.transpose([rel_charges] * 3), axis=0)
                com = _get_weighted_center(poss, [atom.mass for atom in atoms])
                idx = int((com[2] - edges[0]) / dz)
                name_z_dipoles_dict[name][idx].append(dipole)

    for name, z_dipoles in name_z_dipoles_dict.items():
        name_z_dipole_dict[name] = np.array([np.mean(dipoles, axis=0)
                                             if len(dipoles) > n_frame else np.array([0, 0, 0])
                                             for dipoles in z_dipoles])

    name_column_dict = {'z': z_array}
    fig, ax = plt.subplots()
    ax.set(xlabel='z (nm)', ylabel='mean dipole (e*nm)')
    for name, z_dipole in name_z_dipole_dict.items():
        ax.plot(z_array, z_dipole[:, 2], label='dipoleZ - ' + name)
        name_column_dict.update({
            'dipX-' + name: z_dipole[:, 0],
            'dipY-' + name: z_dipole[:, 1],
            'dipZ-' + name: z_dipole[:, 2],
        })

    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-dipole.png')

    print_data_to_file(name_column_dict, f'{args.output}-dipole.txt')


def voltage():
    charges = np.zeros(n_bin)
    top_atom_charges = np.array([atom.charge for atom in top.atoms], dtype=np.float32)
    n_frame = 0
    for i in range(args.begin, args.end, args.skip):
        n_frame += 1
        frame = trj.read_frame(i)
        sys.stdout.write('\r    frame %i' % i)

        positions = frame.positions
        if args.reverse:
            positions[:, 2] = elecd_l + elecd_r - positions[:, 2]

        z_to_left = positions[:, 2] - elecd_l
        ids_ils = np.where((z_to_left > 0.1) & (z_to_left < lz - 0.1))[0]
        z_ils = positions[ids_ils, 2]
        q_ils = frame.charges[ids_ils] if frame.has_charge else top_atom_charges[ids_ils]
        i_bin_ils = ((z_ils - edges[0]) / dz).astype(int)
        for i, i_bin in enumerate(i_bin_ils):
            charges[i_bin] += q_ils[i]

    charges /= n_frame
    ### add back charges on electrodes
    if args.charge != 0:
        q_electrode = args.charge * MILLI * area * NANO * NANO / ELEMENTARY_CHARGE
        idx_cat = int((cathode - edges[0]) / dz)
        idx_ano = int((anode - edges[0]) / dz)
        charges[idx_cat] += q_electrode
        charges[idx_ano] -= q_electrode

    charges_cumulative = np.cumsum(charges)
    charges /= area * dz  # e/nm^3

    e_field = itg.cumtrapz(charges, dx=dz, initial=0) \
              * ELEMENTARY_CHARGE / NANO ** 2 / VACUUM_PERMITTIVITY
    voltage = -itg.cumtrapz(e_field, dx=dz, initial=0) * NANO

    name_column_dict = {'z'    : z_array,
                        'rho_q': charges,
                        'cum_q': charges_cumulative,
                        'EF'   : e_field,
                        'V'    : voltage}
    print_data_to_file(name_column_dict, f'{args.output}-voltage.txt')

    fig, ax = plt.subplots()
    ax.set(xlabel='z (nm)', ylabel='charge density (e/nm$^3$)')
    ax.plot(z_array, charges)
    ax.plot(z_array, [0] * n_bin, '--')
    fig.tight_layout()
    fig.savefig(f'{args.output}-charge.png')

    fig, ax = plt.subplots()
    ax.set(xlabel='z (nm)', ylabel='cumulative charges (e)')
    ax.plot(z_array, charges_cumulative)
    ax.plot(z_array, [0] * n_bin, '--')
    fig.tight_layout()
    fig.savefig(f'{args.output}-charge_cumulative.png')

    fig, ax = plt.subplots()
    ax.set(xlabel='z (nm)', ylabel='electric field (V/m)')
    ax.plot(z_array, e_field)
    ax.plot(z_array, [0] * n_bin, '--')
    fig.tight_layout()
    fig.savefig(f'{args.output}-efield.png')

    fig, ax = plt.subplots()
    ax.set(xlabel='z (nm)', ylabel='voltage (V)')
    ax.plot(z_array, voltage)
    ax.plot(z_array, [0] * n_bin, '--')
    fig.tight_layout()
    fig.savefig(f'{args.output}-voltage.png')


def charge_petersen():
    _conv = ELEMENTARY_CHARGE / area / NANO ** 2 * 1000  # convert from charge (e) to charge density (mC/m^2)
    top_atom_charges = np.array([atom.charge for atom in top.atoms], dtype=np.float32)
    frame_list = []
    q_left_list = []
    for i in range(args.begin, args.end, args.skip):
        frame_list.append(i)
        frame = trj.read_frame(i)
        sys.stdout.write('\r    frame %i' % i)

        positions = frame.positions
        if args.reverse:
            positions[:, 2] = elecd_l + elecd_r - positions[:, 2]

        z_to_left = positions[:, 2] - elecd_l
        ids_ils = np.where((z_to_left > 0.1) & (z_to_left < lz - 0.1))[0]
        z_to_left_ils = z_to_left[ids_ils]
        q_ils = frame.charges[ids_ils] if frame.has_charge else top_atom_charges[ids_ils]
        q_left = np.sum(q_ils * z_to_left_ils) / lz + \
                 args.voltage * area / lz * VACUUM_PERMITTIVITY / ELEMENTARY_CHARGE * NANO
        q_left_list.append(q_left)

    charge_densities = np.array(q_left_list) * _conv
    voltages_surface = charge_densities / 1000 * lz * NANO / VACUUM_PERMITTIVITY

    print('\n%-8s %10s %10s %10s %10s' % ('n_frame', 'rho_q', 'var_rho_q', 'q_atom', 'V_surface'))
    print('%-8i %10.4f %10.4f %10.6f %10.4f' % (len(q_left_list), charge_densities.mean(),
                                                charge_densities.var(),
                                                charge_densities.mean() / _conv / len(ids_cathode),
                                                voltages_surface.mean()))
    name_column_dict = {'frame'    : frame_list,
                        'rho_q'    : charge_densities,
                        'V_surface': voltages_surface}
    print_data_to_file(name_column_dict, f'{args.output}-surface_charge.txt')


def charge_fluctuated():
    _conv = ELEMENTARY_CHARGE / area / NANO ** 2 * 1000  # convert from charge (e) to charge density (mC/m^2)
    qtot_list = []
    print('%-8s %10s %10s' % ('step', 'density_q', 'q_atom'))
    for i in range(args.begin, args.end, args.skip):
        frame = trj.read_frame(i)
        if not frame.has_charge:
            raise Exception('Fluctuated charge analysis requires charge information in trajectory')
        qtot = sum(frame.charges[ids_cathode])
        qtot_list.append(qtot)
        print('%-8i %10.4f %10.6f' % (i, qtot * _conv, qtot / len(ids_cathode)))

    print('\n%-8s %10s %10s %10s' % ('n_step', 'rho_q', 'q_atom', 'V_surface'))
    print('%-8i %10.4f %10.6f %10.4f' % (
        len(qtot_list), np.mean(qtot_list) * _conv, np.mean(qtot_list) / len(ids_cathode),
        np.mean(qtot_list) * _conv / 1000 * lz * NANO / VACUUM_PERMITTIVITY))


def permittivity():
    '''
    Calculate the static relative permittivity from the fluctuation of dipole
    \eps_r = 1 + (<M^2> - <M>^2) / (3VkT 4\pi\eps_0)
    '''
    top_atom_charges = np.array([atom.charge for atom in top.atoms], dtype=np.float32)
    frame_list = []
    dipoles_scalar: [float] = []
    for i in range(args.begin, args.end, args.skip):
        frame = trj.read_frame(i)
        sys.stdout.write('\r    frame %i' % i)

        frame_list.append(i)
        charges = frame.charges if frame.has_charge else top_atom_charges
        dipole = np.sum(frame.positions * np.transpose([charges] * 3), axis=0)
        dipoles_scalar.append(np.sqrt(dipole.dot(dipole)))

    eps_list = [1 + np.var(dipoles_scalar[: i + 1]) * (ELEMENTARY_CHARGE * NANO) ** 2 \
                / (3 * frame0.cell.volume * NANO ** 3) \
                / (8.314 * 298 / AVOGADRO) \
                / (4 * PI * VACUUM_PERMITTIVITY)
                for i in range(len(frame_list))]

    print('\nRelative permittivity = ', eps_list[-1])

    name_column_dict = {'frame'       : frame_list,
                        'dipole'      : dipoles_scalar,
                        'permittivity': eps_list}
    print_data_to_file(name_column_dict, f'{args.output}-permittivity.txt')


if __name__ == '__main__':
    if args.cmd == 'dist':
        distribution()
    elif args.cmd == 'diffuse':
        diffusion()
    elif args.cmd == 'voltage':
        voltage()
    elif args.cmd == 'drude':
        drude_dipole()
    elif args.cmd == 'dipole':
        dipole()
    elif args.cmd == 'charge':
        charge_petersen()
    elif args.cmd == 'fluctq':
        charge_fluctuated()
    elif args.cmd == 'permit':
        permittivity()
