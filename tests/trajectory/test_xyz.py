#!/usr/bin/env python3

import pytest
from mstools.trajectory import Trajectory
from mstools.topology import Topology

import os

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    pass


def test_write():
    top = Topology.open(cwd + '/files/100-SPCE.psf')
    xtc = Trajectory.open(cwd + '/files/100-SPCE.xtc')
    xyz = Trajectory.open(cwd + '/files/xtc-out.xyz', 'w')

    for i in range(xtc.n_frame):
        frame = xtc.read_frame(i)
        xyz.write_frame(frame, top, subset=list(range(150, 300)))
