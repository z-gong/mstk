#!/usr/bin/env python3

import os
import pytest
from mstools.topology import Topology
from mstools.trajectory import Trajectory

cwd = os.path.dirname(os.path.abspath(__file__))

top = Topology.open(cwd + '/files/100-SPCE.psf')
omm_top = top.to_omm_topology()
top = Topology()
top.init_from_omm_topology(omm_top)


def test_read():
    xtc = Trajectory.open(cwd + '/files/100-SPCE.xtc')
    assert xtc.n_atom == 300
    assert xtc.n_frame == 4

    frame = xtc.read_frame(0)
    assert frame.has_velocity == False
    assert frame.has_charge == False
    assert pytest.approx(frame.positions[0], abs=1E-4) == [1.283, 1.791, 1.658]
    assert pytest.approx(frame.positions[-1], abs=1E-4) == [2.726, 2.750, 1.876]

    frame2, frame1 = xtc.read_frames([2, 1])
    assert pytest.approx(frame1.positions[0], abs=1E-4) == [1.335, 1.775, 1.818]
    assert pytest.approx(frame2.positions[-1], abs=1E-4) == [2.529, 0.136, 1.780]


def test_write():
    gro = Trajectory.open(cwd + '/files/100-SPCE.gro')
    xtc = Trajectory.open(cwd + '/files/gro-out.xtc', 'w')

    for i in range(gro.n_frame):
        frame = gro.read_frame(i)
        xtc.write_frame(frame)
