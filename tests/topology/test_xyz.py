#!/usr/bin/env python3

import os
import tempfile
import filecmp
import pytest
import shutil
from mstk.topology import Topology, Molecule

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    xyz = Topology.open(cwd + '/files/urea.xyz')
    assert xyz.n_atom == 8
    assert xyz.cell.volume == 0

    mol = xyz.molecules[0]
    assert mol.name == 'urea'

    atom = xyz.atoms[-1]
    assert atom.name == 'H8'
    assert atom.type == 'HU'
    assert atom.symbol == 'H'
    assert pytest.approx(atom.position, abs=1E-6) == [-0.194866, -0.115100, -0.158814]
    assert atom.has_position


def test_write():
    tmpdir = tempfile.mkdtemp()
    zmat = Topology.open(cwd + '/files/Im11.zmat')
    tmp = os.path.join(tmpdir, 'zmat-out.xyz')
    zmat.write(tmp)
    assert filecmp.cmp(tmp, cwd + '/files/baselines/zmat-out.xyz')
