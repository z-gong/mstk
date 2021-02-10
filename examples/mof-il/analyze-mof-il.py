#!/usr/bin/env python3

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mstools.topology import Topology
from mstools.trajectory import Trajectory

top = Topology.open('topol.psf')

idx_zn1 = np.array([atom.id for atom in top.atoms if atom.symbol == 'Zn' and atom.molecule.id == 0])
idx_zn2 = np.array([atom.id for atom in top.atoms if atom.symbol == 'Zn' and atom.molecule.id == 1])
idx_p = np.array([atom.id for atom in top.atoms if atom.symbol == 'P'])
idx_n = np.array([atom.id for atom in top.atoms if atom.symbol == 'N' and atom.molecule.name == 'NTF2'])
idx_n = np.array([atom.id for atom in top.atoms])

print(len(idx_zn1), len(idx_zn2), len(idx_p), len(idx_n))

trj = Trajectory.open('il-3.65.gro')


def calc_rsq(xyz1, xyz2):
    return np.sum((xyz2 - xyz1) ** 2)


def mesh_density():
    pass


dr = 0.1
max_r = 3
max_z = 5
min_z = -5
mesh = np.zeros([int(max_r / dr), int((max_z - min_z) / dr)])
print(mesh)

r_list = []
z_list = []
for i in range(trj.n_frame):
    sys.stdout.write('\r\t %i / %i' % (i + 1, trj.n_frame))

    frame = trj.read_frame(i)

    rho_p = len(idx_p) / frame.cell.volume
    rho_n = len(idx_n) / frame.cell.volume

    pos_zn1 = frame.positions[idx_zn1]
    com_zn1 = np.mean(pos_zn1, axis=0)

    pos_zn2 = frame.positions[idx_zn2]
    com_zn2 = np.mean(pos_zn2, axis=0)

    com = np.mean([com_zn1, com_zn2], axis=0)
    vec = com_zn2 - com
    vec_unit = vec / np.sqrt(vec.dot(vec))

    print(com_zn1, com_zn2, vec, vec_unit)

    for idx in idx_n:
        pos = frame.positions[idx]
        v = pos - com
        _z = v.dot(vec_unit)
        _r = np.sqrt(v.dot(v) - _z ** 2)

        r_list.append(_r)
        z_list.append(_z)

        if _r < max_r and _z < max_z and _z > min_z:
            mesh[int(_r / dr)][int((_z - min_z) / dr)] += 1 / (2 * np.pi * _r * dr * dr) / rho_n

print(mesh / trj.n_frame)
quit()

fig, ax = plt.subplots()
ax.set(xlim=[0, 3], ylim=[-6, 6])
ax.scatter(r_list, z_list, s=[5.0] * len(r_list), alpha=0.2)
plt.show()
