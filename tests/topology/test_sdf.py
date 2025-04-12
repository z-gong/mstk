#!/usr/bin/env python3

import os
import pytest
import tempfile
import filecmp
from mstk.topology import Topology
from mstk.chem.constant import RAD2DEG

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    top = Topology.open(cwd + '/files/CH3NH3+.sdf')
    assert top.n_molecule == 1

    mol = top.molecules[0]
    assert mol.name == 'CH3NH3+'
    assert mol.n_atom == 8

    atom = top.atoms[1]
    assert atom.name == 'N2'
    assert atom.type == 'UNK'
    assert atom.symbol == 'N'
    assert pytest.approx(atom.position, abs=1E-6) == [0.07355, 0.00131, -0.00126]
    assert atom.has_position
    assert atom.formal_charge == 1


def test_read_V3000():
    top = Topology.open(cwd + '/files/Structures.sdf')
    assert top.n_atom == 473
    assert top.n_molecule == 3
    assert top.n_residue == 3
    assert top.n_bond == 476
    assert top.n_improper == 80
    assert pytest.approx(top.cell.lengths, abs=1E-6) == [0, 0, 0]
    assert pytest.approx(top.cell.angles * RAD2DEG, abs=1E-6) == [90, 90, 90]

    mol = top.molecules[0]
    assert mol.name == 'C4OH'
    assert mol.n_atom == 15
    assert mol.n_residue == 1

    atom = mol.atoms[0]
    assert atom.name == 'O1'
    assert atom.type == 'UNK'
    assert atom.symbol == 'O'
    assert atom.charge == 0.0
    assert pytest.approx(atom.position, abs=1E-6) == [-0.0460000, 0.1301076, 0.0000000]
    assert atom.has_position
    assert atom.molecule.name == 'C4OH'
    assert atom.residue.name == 'C4OH'

    mol = top.molecules[1]
    atom = mol.atoms[-4]
    assert atom.name == 'O5'
    assert atom.type == 'UNK'
    assert atom.symbol == 'O'
    assert atom.formal_charge == -1.0
    assert atom.charge == 0.0
    assert atom.molecule.name == 'CSO3-'
    assert atom.residue.name == 'CSO3-'

    mol = top.molecules[2]
    assert mol.name == 'amorphous poly(AA-HDO)'
    assert mol.n_atom == 90 * 5
    assert mol.n_bond == 91 * 5
    assert mol.n_residue == 1

    atom = mol.atoms[89]
    assert atom.name == 'H90'
    assert atom.type == 'UNK'
    assert atom.symbol == 'H'
    assert atom.charge == 0.0
    assert pytest.approx(atom.position, abs=1E-6) == [0.6584924, 1.3906516, 1.9601327]
    assert atom.has_position
    assert atom.molecule.name == 'amorphous poly(AA-HDO)'
    assert atom.residue.name == 'amorphous_poly(AA-HDO)'
