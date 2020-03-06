#!/usr/bin/env python3

from mstools.topology import Psf

psf = Psf('test-psf.psf')
assert psf.n_atom == 27492
assert psf.n_molecule == 1004
assert psf.is_drude
assert psf.n_bond == 29488
assert psf.n_angle == 35480
assert psf.n_dihedral == 26000
assert psf.n_improper == 0
assert psf.bonds[-1].atom1.name == 'S11'
assert psf.angles[-1].atom2.name == 'S11'
assert psf.dihedrals[-1].atom3.name == 'N5'

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

pout = Psf('out.psf', 'w')
pout.is_drude = psf.is_drude
pout.init_from_molecules(psf.molecules)
pout.write()
pout.close()
