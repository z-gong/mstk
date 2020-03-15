#!/usr/bin/env python3

import numpy as np
from mstools.topology import Topology

import os
cwd = os.path.dirname(os.path.abspath(__file__))
zmat = Topology.open(cwd + '/files/im11.zmat')
print(zmat)

def test_zmat():
    assert zmat.n_atom == 16
    assert zmat.n_bond == 16
    assert zmat.n_angle == 27
    assert zmat.n_dihedral == 32
    assert zmat.n_improper == 5

    atom = zmat.atoms[0]
    assert all(atom.position - np.array([0., 0., 0.]) < 1E-5)
    assert atom.name == 'N1'
    assert atom.type == 'NA'
    atom = zmat.atoms[-1]
    assert all(atom.position - np.array([0.33690, 0.22495, -0.08898]) < 1E-5)
    assert atom.name == 'H16'
    assert atom.type == 'H1'
