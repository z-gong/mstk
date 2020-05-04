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
    dcd = Trajectory.open(cwd + '/files/100-SPCE.dcd')
    assert dcd.n_atom == 300
    assert dcd.n_frame == 3


def test_write():
    gro = Trajectory.open(cwd + '/files/100-SPCE.gro')
    dcd = Trajectory.open(cwd + '/files/gro-out.dcd', 'w')

    for i in range(gro.n_frame):
        frame = gro.read_frame(i)
        dcd.write_frame(frame)
