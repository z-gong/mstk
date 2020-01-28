#!/usr/bin/env python3

import sys
import argparse
import math
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

from mstools.utils import histogram
from mstools.trajectory.topology import Atom, Molecule, Topology
from mstools.trajectory.lammps import LammpsData, LammpsTrj

parser = argparse.ArgumentParser()
parser.add_argument('cmd', choices=['dist', 'diffuse', 'voltage'], help='The property to analyze')
parser.add_argument('-d', '--data', required=True, type=str, help='Lammps data file for topology information')
parser.add_argument('-i', '--input', required=True, type=str,
                    help='Lammps dump trajectory file for atomic positions and charges')
parser.add_argument('-o', '--output', required=True, type=str, help='Output prefix')
parser.add_argument('-b', '--begin', default=0, type=int, help='first frame to output')
parser.add_argument('-e', '--end', default=-1, type=int, help='last frame to output')
parser.add_argument('--skip', default=1, type=int, help='skip frames between output')
args = parser.parse_args()

top = LammpsData(args.data)
trj = LammpsTrj(args.input)

print('Topology info: ', top.n_atom, 'atoms;', top.n_molecule, 'molecules')
print('Trajectory info: ', trj.n_atom, 'atoms;', trj.n_frame, 'frames')

if (top.n_atom != trj.n_atom):
    raise Exception('Number of atoms in topology and trajectory files do not match')

if args.end > trj.n_frame or args.end == -1:
    args.end = trj.n_frame


def _get_atoms(mol: Molecule, names: [str]):
    return [list(filter(lambda x: x.name == name, mol.atoms))[0] for name in names]


def _get_com_position(positions, atoms):
    atom_masses = np.array([atom.mass for atom in atoms])
    com_mass = sum(atom_masses)
    atom_ids = [atom.id for atom in atoms]
    atom_positions = positions[atom_ids]
    com_position = np.sum(atom_positions * np.transpose(np.array([atom_masses] * 3)), axis=0) / com_mass
    return com_position


def distribution():
    bins = np.linspace(0, 70, 701)
    dz = bins[1] - bins[0]
    z_array = (bins[1:] + bins[:-1]) / 2
    charge_array = np.array([0.] * len(z_array))
    z_ring_list = []
    z_ring_com_list = []
    z_tail_list = []
    z_dca_list = []
    z_dca_com_list = []
    theta_ring_list = []
    theta_tail_list = []
    theta_dca_list = []

    frame = trj.read_frame(0)
    area = frame.box[0] * frame.box[1]

    n_frame = 0
    for i in range(args.begin, args.end, args.skip):
        n_frame += 1
        frame = trj.read_frame(i)
        sys.stdout.write('\r    step %i' % frame.step)
        positions = frame.positions
        for mol in top.molecules:
            if mol.name == 'c8c1im+':
                ring_atoms = _get_atoms(mol, ['N1', 'C3', 'N5', 'C7', 'C9', 'C11', 'C14'])
                z_ring_list += [positions[atom.id][2] for atom in ring_atoms]
                ring_com = _get_com_position(positions, ring_atoms)
                z_ring_com_list.append(ring_com[2])
                tail_atoms = _get_atoms(mol, ['C21', 'C25', 'C29', 'C33', 'C37', 'C41', 'C45'])
                z_tail_list += [positions[atom.id][2] for atom in tail_atoms]
                if ring_com[2] > 64.5:
                    a1, a2 = _get_atoms(mol, ['C21', 'C45'])
                    vec_tail = positions[a2.id] - positions[a1.id]
                    vec_proj = np.array([vec_tail[0], vec_tail[1], 0])
                    theta = np.arccos(vec_tail.dot(vec_proj) / np.sqrt(vec_tail.dot(vec_tail)) / np.sqrt(
                        vec_proj.dot(vec_proj))) * 180 / np.pi
                    theta_tail_list.append(theta)

                    a1, a2, a3 = _get_atoms(mol, ['C3', 'C7', 'C9'])
                    vec_12 = positions[a2.id] - positions[a1.id]
                    vec_13 = positions[a3.id] - positions[a1.id]
                    vec_normal = np.cross(vec_12, vec_13)
                    vec_proj = np.array([vec_normal[0], vec_normal[1], 0])
                    theta = np.arccos(vec_normal.dot(vec_proj) / np.sqrt(vec_normal.dot(vec_normal)) / np.sqrt(
                        vec_proj.dot(vec_proj))) * 180 / np.pi
                    theta_ring_list.append(90 - theta)
            if mol.name == 'dca-':
                dca_atoms = _get_atoms(mol, ['N1', 'C3', 'N5', 'C7', 'N9'])
                z_dca_list += [positions[atom.id][2] for atom in dca_atoms]
                dca_com = _get_com_position(positions, dca_atoms)
                z_dca_com_list.append(dca_com[2])
                if dca_com[2] < 5:
                    a1, a2 = _get_atoms(mol, ['N5', 'N9'])
                    vec_dca = positions[a2.id] - positions[a1.id]
                    vec_proj = np.array([vec_dca[0], vec_dca[1], 0])
                    theta = np.arccos(vec_dca.dot(vec_proj) / np.sqrt(vec_dca.dot(vec_dca)) / np.sqrt(
                        vec_proj.dot(vec_proj))) * 180 / np.pi
                    theta_dca_list.append(theta)

        for atom in top.atoms:
            if atom.molecule.name not in ['c8c1im+', 'dca-']:
                continue
            z = frame.positions[atom.id][2]
            q = frame.charges[atom.id]
            z_idx = math.floor((z - bins[0]) / dz)
            charge_array[z_idx] += q

    print('')

    fig, ax = plt.subplots()
    ax.set(xlim=[0, 70], xlabel='z (A)', ylabel='particle density (/$A^3$)')
    x, y = histogram(z_ring_list, bins=bins)
    ax.plot(x, y / area / dz / n_frame, label='c8c1im+ ring')
    x, y = histogram(z_tail_list, bins=bins)
    ax.plot(x, y / area / dz / n_frame, label='c8c1im+ tail')
    x, y = histogram(z_dca_list, bins=bins)
    ax.plot(x, y / area / dz / n_frame, label='dca-')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-dist.png')
    plt.cla()

    ax.set(xlim=[0, 70], xlabel='z (A)', ylabel='particle density (/$A^3$)')
    x, y = histogram(z_ring_com_list, bins=bins)
    ax.plot(x, y / area / dz / n_frame, label='c8c1im+ ring COM')
    x, y = histogram(z_dca_com_list, bins=bins)
    ax.plot(x, y / area / dz / n_frame, label='dca- COM')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-dist-com.png')
    plt.cla()

    ax.set(xlim=[0, 70], xlabel='z (A)', ylabel='charge density (e/$A^3$)')
    ax.plot(z_array, charge_array / area / dz / n_frame, label='c8c1im+ dca-')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-charge.png')
    plt.cla()

    ax.set(xlim=[0, 90], ylim=[0, 500], xlabel='theta (A)', ylabel='probability')
    x, y = histogram(theta_ring_list, bins=np.linspace(0, 90, 91))
    ax.plot(x, y, label='c8c1im+ ring direction')
    x, y = histogram(theta_tail_list, bins=np.linspace(0, 90, 91))
    ax.plot(x, y, label='c8c1im+ tail direction')
    x, y = histogram(theta_dca_list, bins=np.linspace(0, 90, 91))
    ax.plot(x, y, label='dca- end-end direction')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{args.output}-angle.png')
    plt.cla()


def diffusion():
    ring_atoms_list = []
    dca_atoms_list = []

    frame = trj.read_frame(args.begin)
    positions = frame.positions
    for mol in top.molecules:
        if mol.name == 'c8c1im+':
            ring_atoms = _get_atoms(mol, ['N1', 'C3', 'N5', 'C7', 'C9', 'C11', 'C14'])
            ring_com = _get_com_position(positions, ring_atoms)
            if ring_com[2] > 54.5:
                ring_atoms_list.append(ring_atoms)

        if mol.name == 'dca-':
            dca_atoms = _get_atoms(mol, ['N1', 'C3', 'N5', 'C7', 'N9'])
            dca_com = _get_com_position(positions, dca_atoms)
            if dca_com[2] < 15:
                dca_atoms_list.append(dca_atoms)

    t_list = []
    ring_com_z_list = [[] for i in ring_atoms_list]
    dca_com_z_list = [[] for i in dca_atoms_list]

    for i in range(args.begin, args.end, args.skip):
        frame = trj.read_frame(i)
        sys.stdout.write('\r    step %i' % frame.step)

        t_list.append(frame.step / 1000)  # ps
        for i, ring_atoms in enumerate(ring_atoms_list):
            ring_com = _get_com_position(frame.positions, ring_atoms)
            ring_com_z_list[i].append(ring_com[2])
        for i, dca_atoms in enumerate(dca_atoms_list):
            dca_com = _get_com_position(frame.positions, dca_atoms)
            dca_com_z_list[i].append(dca_com[2])

    print('')

    fig, ax = plt.subplots()
    ax.set(ylim=[54.5, 69.5], xlabel='time (ps)', ylabel='z (A)')
    for z_list in ring_com_z_list:
        ax.plot(t_list, z_list)
    fig.tight_layout()
    fig.savefig(f'{args.output}-diffusion-ring.png')
    plt.cla()

    ax.set(ylim=[0, 15], xlabel='time (ps)', ylabel='z (A)')
    for z_list in dca_com_z_list:
        ax.plot(t_list, z_list)
    fig.tight_layout()
    fig.savefig(f'{args.output}-diffusion-dca.png')
    plt.cla()


def voltage():
    frame = trj.read_frame(0)
    area = frame.box[1] * frame.box[2]
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
            z = frame.positions[k][2]
            i_bin = math.floor((z - frame.zlo) / dz)
            i_bin = min(i_bin, n_bin - 1)
            charge = frame.charges[k] if frame.has_charge else atom.charge
            charges[i_bin] += charge

    charges_cumulative = np.cumsum(charges) / n_frame
    charges /= area * dz * n_frame  # e/A^3

    eps0 = 8.854E-12
    q0 = 1.6E-19
    Ang = 1E-10

    for i in range(1, n_bin):
        s = 0
        for j in range(0, i + 1):
            s += dz * (dz * (i - j)) * charges[j]
        voltage[i] = -s / eps0 * q0 / Ang

    fig = plt.figure(figsize=[6.4, 12.0])
    ax1 = fig.add_subplot('311')
    ax1.set(ylabel='charge density (e/A$^3$)')
    ax1.plot(bins, charges)
    ax1.plot(bins, [0] * n_bin, '--')

    ax2 = fig.add_subplot('312')
    ax2.set(ylabel='cumulative charges (e)')
    ax2.plot(bins, charges_cumulative)
    ax2.plot(bins, [0] * n_bin, '--')

    ax3 = fig.add_subplot('313')
    ax3.set(ylabel='voltage (V)')
    ax3.plot(bins, voltage)
    ax3.plot(bins, [0] * n_bin, '--')
    fig.tight_layout()
    fig.savefig(f'{args.output}-voltage.png')


if __name__ == '__main__':
    if args.cmd == 'dist':
        distribution()
    elif args.cmd == 'diffuse':
        diffusion()
    elif args.cmd == 'voltage':
        voltage()
