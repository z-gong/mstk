#!/usr/bin/env python3

from mstools.trajectory import Gro

import os
cwd = os.path.dirname(os.path.abspath(__file__))
gro = Gro(cwd + '/files/dump.gro')

def test_gro():
    assert gro.n_atom == 1000
    assert gro.n_frame == 3

    frame = gro.read_frame(0)
    assert frame.has_velocity == True
    assert all(frame.positions[0] == [1.971, 3.328, 4.796])

    frame = gro.read_frame(2)
    assert all(frame.box == [4.91065, 4.91067, 4.91066])
    assert all(frame.velocities[-1] == [0.1388, -0.2205,  0.3017])
