#!/usr/bin/env python3

import math
import sys
import numpy as np

from mstools.trajectory.lammpstrj import LammpsTrj
from mstools.trajectory.pdb import PDB

trj = LammpsTrj('see.lammpstrj')
print(trj.n_atom, 'atoms;', trj.n_frame, 'frames')

pdb = PDB('dump.pdb', 'w')
for i in range(10, trj.n_frame,5):
    sys.stdout.write('\r    %i' % i)
    frame = trj.read_frame(i)
    pdb.write_frame(frame)
pdb.close()
