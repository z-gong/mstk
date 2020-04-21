#!/usr/bin/env python3

import pytest
from mstools.trajectory import Gro

import os
cwd = os.path.dirname(os.path.abspath(__file__))

def test_read():
    gro = Gro(cwd + '/files/1000-W.gro')
    assert gro.n_atom == 1000
    assert gro.n_frame == 3

    frame = gro.read_frame(0)
    assert frame.has_velocity == True
    assert pytest.approx(frame.positions[0], abs=1E-6) == [1.971, 3.328, 4.796]

    frame = gro.read_frame(2)
    assert frame.cell.is_rectangular
    assert pytest.approx(frame.cell.size, abs=1E-6) == [4.91065, 4.91067, 4.91066]
    assert pytest.approx(frame.velocities[-1], abs=1E-6) == [0.1388, -0.2205,  0.3017]
