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


def test_write():
    xtc = Trajectory.open(cwd + '/files/100-SPCE.xtc')
    gro = Trajectory.open(cwd + '/files/xtc-out.gro', 'w')
    dcd = Trajectory.open(cwd + '/files/xtc-out.dcd', 'w')
    out = Trajectory.open(cwd + '/files/xtc-out.xtc', 'w')

    frame = xtc.read_frame(0)
    assert frame.has_velocity == False
    gro.write_frame(frame, top, subset=list(range(100)))
    dcd.write_frame(frame)
    out.write_frame(frame)

    frame = xtc.read_frame(2)
    assert frame.cell.is_rectangular
    gro.write_frame(frame, top, subset=list(range(100)))
    dcd.write_frame(frame)
    out.write_frame(frame)
