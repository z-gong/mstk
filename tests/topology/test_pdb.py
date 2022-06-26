#!/usr/bin/env python3

import os
import pytest
import tempfile
import filecmp
import shutil
from mstk.topology import Topology, Molecule

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    pdb = Topology.open(cwd + '/files/test.pdb')
    assert pdb.n_atom == 16
    assert pytest.approx(pdb.cell.volume, abs=1E-6) == 0.06

    mol = pdb.molecules[0]
    assert mol.name == 'c1c1'

    atom = pdb.atoms[-1]
    assert atom.name == 'H16'
    assert atom.type == ''
    assert atom.symbol == 'H'
    assert pytest.approx(atom.position, abs=1E-6) == [0.3369, 0.2249, -0.0890]
    assert atom.has_position


def test_write():
    tmpdir = tempfile.mkdtemp()
    zmat = Topology.open(cwd + '/files/Im11.zmat')
    tmp = os.path.join(tmpdir, 'zmat-out.pdb')
    zmat.write(tmp)
    assert filecmp.cmp(tmp, cwd + '/files/baselines/zmat-out.pdb')
    shutil.rmtree(tmpdir)
