#!/usr/bin/env python3

from mstools.topology import Psf

psf = Psf('test-psf.psf')
assert psf.n_atom == 27492
assert psf.n_molecule == 1004
assert psf.is_drude

atom = psf.atoms[0]
assert atom.id == 0
assert atom.molecule.id == 0
assert atom.type == 'MoS'
assert atom.name == 'Mo1'
assert atom.charge == 4.924306
assert atom.mass == 95.5370

atom = psf.atoms[-1]
assert atom.id == 27491
assert atom.molecule.id == 1003
assert atom.type == 'D'
assert atom.name == 'D18'
assert atom.charge == -1.371919
assert atom.mass == 0.4

mol = psf.molecules[0]
assert mol.id == 0
assert mol.name == 'MoS'

mol = psf.molecules[-1]
assert mol.id == 1003
assert mol.name == 'fsi'
