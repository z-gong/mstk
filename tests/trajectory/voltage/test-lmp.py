#!/usr/bin/env python3

import math
import sys
import numpy as np

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

from mstools.trajectory.lammpstrj import LammpsTrj
from mstools.molio import Atom

trj = LammpsTrj('dump.lammpstrj')
print(trj.n_atom, 'atoms;', trj.n_frame, 'frames')

frame = trj.read_frame(0)
area = frame.box[1] * frame.box[2]
dz = 0.1
n_bin = math.ceil(frame.box[2] / dz)
bins = np.array([frame.zlo + dz * (i + 0.5) for i in range(n_bin)])
charges = np.array([0.] * n_bin)
voltage = np.array([0.] * n_bin)

n_frame = 0
for i in range(0, trj.n_frame, 2):
    n_frame += 1
    frame = trj.read_frame(i)
    sys.stdout.write('\r    step %i' % frame.step)

    atom: Atom
    for atom in frame.atoms:
        z = atom.position[2]
        i_bin = math.floor((z - frame.zlo) / dz)
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
fig.savefig('test.png')
fig.show()
