#!/usr/bin/env python3

import os
import tempfile
import filecmp
import pytest
import shutil
from mstools.topology import Topology
from mstools.trajectory import Trajectory

cwd = os.path.dirname(os.path.abspath(__file__))


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
    tmpdir = tempfile.mkdtemp()
    gro = Trajectory.open(cwd + '/files/100-SPCE.gro')
    tmp = os.path.join(tmpdir, 'gro-out.xtc')
    xtc = Trajectory.open(tmp, 'w')
    for i in range(gro.n_frame):
        frame = gro.read_frame(i)
        xtc.write_frame(frame)
    xtc.close()
    assert filecmp.cmp(tmp, cwd + '/files/baselines/gro-out.xtc')
    shutil.rmtree(tmpdir)
