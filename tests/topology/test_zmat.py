#!/usr/bin/env python3

import os
import pytest
from mstk.topology import Topology

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    zmat = Topology.open(cwd + '/files/Im11.zmat')
    assert zmat.n_atom == 16
    assert zmat.n_bond == 16
    assert zmat.n_angle == 27
    assert zmat.n_dihedral == 32
    assert zmat.n_improper == 5
    assert zmat.cell.volume == 0

    atom = zmat.atoms[0]
    assert all(atom.position == [0., 0., 0.])
    assert atom.name == 'N1'
    assert atom.type == 'NA'
    atom = zmat.atoms[-1]
    assert pytest.approx(atom.position, abs=1E-4) == [0.3369, 0.2250, -0.0890]
    assert atom.name == 'H16'
    assert atom.type == 'H1'
