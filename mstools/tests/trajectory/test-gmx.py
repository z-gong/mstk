#!/usr/bin/env python3

import math
import sys
import numpy as np

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

from mstools.trajectory.lammpsdump import LammpsDump
from mstools.molio import Atom

import MDAnalysis as mda

u = mda.Universe('test-MoS2-nopol.tpr', 'test-MoS2-nopol-gmx.xtc')
trj = u.trajectory
print(trj.n_atoms, len(trj))

frame = trj[0]
area = frame.dimensions[0] * frame.dimensions[1]
dz = 0.1
n_bin = math.ceil(frame.dimensions[2] / dz)
bins = np.array([dz * (i + 0.5) for i in range(n_bin)])
charges = np.array([0.] * n_bin)
voltage = np.array([0.] * n_bin)

n_frame = 0
for frame in trj:
    n_frame += 1
    print(frame.frame)

    for i, pos in enumerate(frame.positions):
        atom = u.atoms[i]
        z = pos[2]
        i_bin = math.floor(z / dz)
        i_bin = min(i_bin, n_bin - 1)
        charges[i_bin] += atom.charge

charges /= area * dz * n_frame  # e/A^3

eps0 = 8.854E-12
q0 = 1.6E-19
Ang = 1E-10

for i in range(1, n_bin):
    s = 0
    for j in range(0, i + 1):
        s += dz * (dz * (i - j)) * charges[j]
    voltage[i] = -s / eps0 * q0 / Ang

fig = plt.figure(figsize=[6.4, 8.0])
ax1 = fig.add_subplot('211')
ax1.set(ylabel='charge density (e/A$^3$)')
ax1.plot(bins, charges)
ax1.plot(bins, [0] * n_bin, '--')
ax2 = fig.add_subplot('212')
ax2.set(ylabel='voltage (V)')
ax2.plot(bins, voltage)
ax2.plot(bins, [0] * n_bin, '--')
fig.tight_layout()
fig.savefig('test-MoS2-nopol-gmx.png')
fig.show()
