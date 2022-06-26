#!/usr/bin/env python3

import os
import tempfile
import filecmp
import pytest
import shutil
from mstk.trajectory import Trajectory

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    dcd = Trajectory.open(cwd + '/files/100-SPCE.dcd')
    assert dcd.n_atom == 300
    assert dcd.n_frame == 3

    frame = dcd.read_frame(0)
    assert frame.has_velocity == False
    assert frame.has_charge == False
    assert pytest.approx(frame.positions[0], abs=1E-6) == [1.180277, 1.421659, 1.885429]

    frame2, frame1 = dcd.read_frames([2, 1])
    assert pytest.approx(frame1.positions[0], abs=1E-6) == [1.645082, 1.233684, 2.165311]
    assert pytest.approx(frame2.positions[-1], abs=1E-6) == [2.545774, 2.661569, 1.707037]
    assert pytest.approx(frame2.cell.lengths, abs=1E-6) == [2.983427, 2.983427, 2.983427]


def test_write():
    tmpdir = tempfile.mkdtemp()

    gro = Trajectory.open(cwd + '/files/100-SPCE.gro')
    tmp = os.path.join(tmpdir, 'gro-out.dcd')
    dcd = Trajectory.open(tmp, 'w')
    for i in range(gro.n_frame):
        frame = gro.read_frame(i)
        dcd.write_frame(frame)
    dcd.close()
    # It's annoying that DCD files generated are different every time. Compare the size instead of the content.
    assert os.path.getsize(tmp) == os.path.getsize(cwd + '/files/baselines/gro-out.dcd')

    shutil.rmtree(tmpdir)
