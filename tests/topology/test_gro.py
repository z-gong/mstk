#!/usr/bin/env python3

import os
import pytest
import tempfile
import filecmp
import shutil
from mstk.topology import Topology, Molecule

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    gro = Topology.open(cwd + '/files/100-SPCE.gro')
    assert gro.n_atom == 300
    assert pytest.approx(gro.cell.volume, abs=1E-6) == 27.0

    assert gro.n_molecule == 100
    mol = gro.molecules[0]
    assert mol.name == 'SPCE'
    assert mol.n_atom == 3

    assert gro.n_residue == 100
    res = gro.residues[-1]
    assert res.name == 'SPCE'
    assert res.n_atom == 3

    atom = gro.atoms[-1]
    assert atom.name == 'H'
    assert atom.type == 'H'
    assert atom.symbol == 'H'
    assert pytest.approx(atom.position, abs=1E-6) == [2.726, 2.750, 1.876]
    assert atom.has_position


def test_read_repeated_resid():
    gro = Topology.open(cwd + '/files/2chain.gro')
    assert gro.n_atom == 26
    assert gro.n_molecule == 6
    assert gro.molecules[0].name == 'INI'
    assert gro.molecules[0].n_atom == 1
    assert gro.molecules[1].name == 'B'
    assert gro.molecules[1].n_atom == 4
    assert gro.molecules[2].name == 'A'
    assert gro.molecules[2].n_atom == 8
    assert gro.molecules[3].name == 'INI'
    assert gro.molecules[3].n_atom == 1
    assert gro.molecules[4].name == 'B'
    assert gro.molecules[4].n_atom == 4
    assert gro.molecules[5].name == 'A'
    assert gro.molecules[5].n_atom == 8


def test_write():
    tmpdir = tempfile.mkdtemp()
    zmat = Topology.open(cwd + '/files/100-SPCE.gro')
    tmp = os.path.join(tmpdir, 'gro-out.gro')
    zmat.write(tmp)
    assert filecmp.cmp(tmp, cwd + '/files/baselines/gro-out.gro')
    shutil.rmtree(tmpdir)
