#!/usr/bin/env python3

import pytest
from mstools.topology import Topology

import os
cwd = os.path.dirname(os.path.abspath(__file__))
xyz = Topology.open(cwd + '/files/urea.xyz')

def test_xyz():
    assert xyz.n_atom == 8
    assert xyz.remark == 'urea'
    atom = xyz.atoms[-1]
    assert atom.name == 'H8'
    assert atom.type == 'HU'
    assert atom.symbol == 'H'
    assert pytest.approx(atom.position, abs=1E-6) == [-0.194866, -0.115100, -0.158814]
    assert atom.has_position
