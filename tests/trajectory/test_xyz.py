#!/usr/bin/env python3

import os
import tempfile
import filecmp
import pytest
import shutil
from mstk.trajectory import Trajectory
from mstk.topology import Topology

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    xyz = Trajectory.open(cwd + '/files/100-SPCE.xyz')
    assert xyz.n_frame == 4
    assert xyz.n_atom == 300
    frame1 = xyz.read_frame(1)
    assert pytest.approx(frame1.positions[150], abs=1E-6) == [1.227000, 1.109000, 2.458000]
    frame2, frame0 = xyz.read_frames([2, 0])
    assert pytest.approx(frame0.positions[150], abs=1E-6) == [1.410000, 1.315000, 2.851000]
    assert pytest.approx(frame2.positions[-1], abs=1E-6) == [2.529000, 0.136000, 1.780000]


def test_write():
    tmpdir = tempfile.mkdtemp()

    top = Topology.open(cwd + '/files/100-SPCE.psf')
    xtc = Trajectory.open(cwd + '/files/100-SPCE.xtc')
    tmp = os.path.join(tmpdir, 'xtc-out.xyz')
    xyz = Trajectory.open(tmp, 'w')
    for i in range(xtc.n_frame):
        frame = xtc.read_frame(i)
        xyz.write_frame(frame, top, subset=list(range(150, 300)))
    xyz.close()
    filecmp.cmp(tmp, cwd + '/files/baselines/xtc-out.xyz')

    shutil.rmtree(tmpdir)
