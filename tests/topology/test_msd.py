#!/usr/bin/env python3

import os
import pytest
import tempfile
import filecmp
from mstools.topology import Topology, Molecule

cwd = os.path.dirname(os.path.abspath(__file__))
tempdir = tempfile.mkdtemp()


def test_read():
    msd = Topology.open(cwd + '/files/im41-sol.msd')
    assert msd.n_atom == 28
    assert msd.n_molecule == 2
    assert msd.n_bond == 27
    assert msd.n_improper == 5
    assert pytest.approx(msd.cell.lengths, abs=1E-6) == [3.0, 2.5, 2.0]
    assert pytest.approx(msd.cell.angles, abs=1E-6) == [90, 90, 60]

    mol = msd.molecules[0]
    assert mol.name == 'UNK'
    assert mol.n_atom == 25

    mol = msd.molecules[1]
    assert mol.name == 'SOL'
    assert mol.n_atom == 3
    assert mol.n_bond == 2

    atom = msd.atoms[0]
    assert atom.name == 'C1'
    assert atom.type == 'C02'
    assert atom.symbol == 'C'
    assert atom.charge == -0.06
    assert pytest.approx(atom.position, abs=1E-6) == [-0.3588655, -0.1294565, -0.0387042]
    assert atom.has_position
    assert atom.molecule.name == 'UNK'

    atom = msd.atoms[-1]
    assert atom.name == 'H28'
    assert atom.type == 'H08'
    assert atom.symbol == 'H'
    assert atom.charge == 0.4238
    assert pytest.approx(atom.position, abs=1E-6) == [0.18287, 0.1195811, 0.0488527]
    assert atom.has_position
    assert atom.molecule.name == 'SOL'
