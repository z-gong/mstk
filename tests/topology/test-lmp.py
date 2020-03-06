#!/usr/bin/env python3

from mstools.topology import LammpsData

lmp = LammpsData('test-lmp.lmp')
assert lmp.n_atom == 20594
assert lmp.n_molecule == 9895
assert lmp.is_drude == False

atom = lmp.atoms[0]
assert atom.id == 0
assert atom.molecule.id == 0
assert atom.molecule.name == 'MoS2'
assert atom.type == 'MoS'
assert atom.name == 'Mo1'
assert atom.charge == 0.0
assert atom.mass == 95.937

atom = lmp.atoms[6879]
assert atom.id == 6879
assert atom.molecule.id == 110
assert atom.molecule.name == 'c8c1im+'
assert atom.type == 'D'
assert atom.name == 'D34'
assert atom.charge == -2.237196
assert atom.mass == 0.4

mol = lmp.molecules[0]
assert mol.id == 0
assert mol.name == 'MoS2'

mol = lmp.molecules[-1]
assert mol.id == 9894
assert mol.name == 'IMG'
