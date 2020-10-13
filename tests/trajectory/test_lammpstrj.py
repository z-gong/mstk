#!/usr/bin/env python3

import pytest
from mstools.trajectory import Trajectory

import os

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    trj = Trajectory.open(cwd + '/files/100-SPCE.lammpstrj')
    assert trj.n_atom == 300
    assert trj.n_frame == 4

    frame = trj.read_frame(0)
    assert frame.has_velocity == False
    assert pytest.approx(frame.positions[0], abs=1E-6) == [2.2545, 1.45688, 1.98961]
    assert frame.has_charge == True
    assert pytest.approx(frame.charges[:4], abs=1E-6) == [0.4238, -0.8476, 0.4238, 0.4238]

    frame = trj.read_frame(1)
    assert frame.cell.is_rectangular
    assert pytest.approx(frame.cell.size, abs=1E-6) == [3.0004316] * 3
    assert pytest.approx(frame.positions[-1], abs=1E-6) == [2.11148, 0.241373, 0.664092]
