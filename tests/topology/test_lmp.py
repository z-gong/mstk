#!/usr/bin/env python3

from mstools.topology import LammpsData

import os
cwd = os.path.dirname(os.path.abspath(__file__))
lmp = LammpsData(cwd + '/files/data.lmp')
print(lmp.box)

def test_lmp():
    assert lmp.n_atom == 75
    assert lmp.n_molecule == 15
    assert lmp.is_drude == False

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

    from mstools.topology import Psf
    pout = Psf(cwd + '/files/lmp-out.psf', 'w')
    pout.is_drude = lmp.is_drude
    pout.init_from_molecules(lmp.molecules)
    pout.write()
    pout.close()
