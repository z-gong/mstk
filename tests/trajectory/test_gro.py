#!/usr/bin/env python3

import tempfile
import filecmp
import pytest
import shutil
from mstools.trajectory import Trajectory
from mstools.topology import Topology

import os

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    gro = Trajectory.open(cwd + '/files/100-SPCE.gro')
    assert gro.n_atom == 300
    assert gro.n_frame == 2

    frame = gro.read_frame(0)
    assert frame.has_velocity == True
    assert pytest.approx(frame.positions[0], abs=1E-6) == [1.283, 1.791, 1.658]

    frame = gro.read_frame(1)
    assert frame.cell.is_rectangular
    assert pytest.approx(frame.cell.size, abs=1E-6) == [3.00906, 3.00906, 3.00906]
    assert pytest.approx(frame.velocities[-1], abs=1E-6) == [-1.0323, 0.5604, -0.3797]


def test_write():
    tmpdir = tempfile.mkdtemp()
    top = Topology.open(cwd + '/files/100-SPCE.psf')
    xtc = Trajectory.open(cwd + '/files/100-SPCE.xtc')

    tmp = os.path.join(tmpdir, 'xtc-out.gro')
    gro = Trajectory.open(tmp, 'w')
    for i in range(xtc.n_frame):
        frame = xtc.read_frame(i)
        gro.write_frame(frame, top, subset=list(range(150, 300)))
    gro.close()
    filecmp.cmp(tmp, cwd + '/files/baselines/xtc-out.gro')

    shutil.rmtree(tmpdir)
