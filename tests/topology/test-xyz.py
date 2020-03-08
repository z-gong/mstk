#!/usr/bin/env python3

import numpy as np
from mstools.topology import Topology

xyz = Topology.open('test-xyz.xyz')

assert xyz.n_atom == 8
assert xyz.remark == 'urea'
atom = xyz.atoms[-1]
assert atom.name == 'H8'
assert atom.type == 'HU'
assert atom.symbol == 'H'
assert all(atom.position - np.array([-0.194866, -0.115100, -0.158144]) < 1E-5)
assert atom.has_position
