#!/usr/bin/env python3

import os
import pytest
from mstk.topology import Topology

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    lmp = Topology.open(cwd + '/files/10-H2O-5-C3H6.lmp')
    assert lmp.n_atom == 75
    assert lmp.n_molecule == 15
    assert lmp.is_drude == False
    assert lmp.cell.is_rectangular
    assert pytest.approx(lmp.cell.lengths, abs=1E-6) == [3.18888, 3.18888, 3.18888]

    atom = lmp.atoms[0]
    assert atom.id == 0
    assert atom.molecule.id == 0
    assert atom.molecule.name == 'WAT'
    assert atom.type == 'o_2w'
    assert atom.name == 'O1'
    assert atom.charge == -0.8476
    assert atom.mass == 15.9994

    atom = lmp.atoms[-1]
    assert atom.id == 74
    assert atom.molecule.id == 14
    assert atom.molecule.name == 'C3H6'
    assert atom.type == 'h_1'
    assert atom.name == 'H9'
    assert atom.charge == 0.06
    assert atom.mass == 1.0079

    mol = lmp.molecules[0]
    assert mol.id == 0
    assert mol.name == 'WAT'

    mol = lmp.molecules[-1]
    assert mol.id == 14
    assert mol.name == 'C3H6'

    assert lmp.has_position
    atom = lmp.atoms[0]
    assert pytest.approx(atom.position, abs=1E-6) == [2.4257, 0.3594, 0.3218]
    atom = lmp.atoms[-1]
    assert pytest.approx(atom.position, abs=1E-6) == [1.6725, 2.1756, 0.5918]


def test_improper():
    lmp = Topology.open(cwd + '/files/10-H2O-5-C3H6.lmp')
    assert lmp.n_improper == 10
    improper = lmp.impropers[0]
    assert improper.atom1.id + 1 == 32
    assert improper.atom2.id + 1 == 34
    assert improper.atom3.id + 1 == 31
    assert improper.atom4.id + 1 == 35

    lmp = Topology.open(cwd + '/files/10-H2O-5-C3H6.lmp', improper_center=3)
    improper = lmp.impropers[0]
    assert improper.atom1.id + 1 == 31
    assert improper.atom2.id + 1 == 32
    assert improper.atom3.id + 1 == 34
    assert improper.atom4.id + 1 == 35
