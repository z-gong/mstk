#!/usr/bin/env python3

import numpy as np
from mstools.trajectory import Gro

import os
cwd = os.path.dirname(os.path.abspath(__file__))
gro = Gro(cwd + '/files/dump.gro')

def test_gro():
    assert gro.n_atom == 1000
    assert gro.n_frame == 3

    frame = gro.read_frame(0)
    assert frame.has_velocity == True
    assert np.all(frame.positions[0] == np.array([1.971, 3.328, 4.796]))

    frame = gro.read_frame(2)
    assert np.all(frame.box == np.array([4.91065, 4.91067, 4.91066]))
    assert np.all(frame.velocities[-1] == np.array([0.1388, -0.2205,  0.3017]))
