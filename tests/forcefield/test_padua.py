#!/usr/bin/env python3

import pytest
from mstools.forcefield import PaduaFFSet

import os

cwd = os.path.dirname(os.path.abspath(__file__))
params = PaduaFFSet(cwd + '/files/oplsaa.ff', cwd + '/files/clp-alpha.ff')


def test_read():
    br = params.atom_types['Br']
    assert br.mass == 79.904
    assert br.charge == -1.0

    cap = params.atom_types['CAP']
    ca = params.atom_types['CA']
    assert cap.eqt_vdw == 'CAP'
    assert cap.eqt_bond == 'CA'
    assert cap.eqt_imp_s == ca.name
    assert ca.eqt_vdw == ca.name
    assert ca.eqt_ang_c == ca.name
    assert ca.eqt_dih_c == 'CA'

    vdw = params.vdw_terms['Br,Br']
    assert vdw.epsilon == 0.37656
    assert vdw.sigma == 4.624 / 10

    bond = params.bond_terms['CT,CT']
    assert bond.length == 1.529 / 10
    assert pytest.approx(bond.k, abs=1E-6) == 2242.6 / 2 * 100

    angle = params.angle_terms['CT,CT,HC']
    assert angle.theta == 110.7
    assert angle.k == 313.8 / 2

    dihedral = params.dihedral_terms['CT,CT,CT,HC']
    assert dihedral.k1 == 0
    assert dihedral.k2 == 0
    assert dihedral.k3 == 1.2552 / 2
    assert dihedral.k4 == 0

    dihedral = params.dihedral_terms['CT,CT,CT,CT']
    assert dihedral.k1 == 5.4392 / 2
    assert dihedral.k2 == -0.2092 / 2
    assert dihedral.k3 == 0.8368 / 2
    assert dihedral.k4 == 0

    improper = params.improper_terms['N,C,CT,CT']
    assert improper.phi == 180
    assert improper.k == 8.368 / 2

    drude = params.polarizable_terms['CT']
    assert drude.mass == 0.4
    assert pytest.approx(drude.k / 100, abs=1E-6) == 4184 / 2
    assert pytest.approx(drude.alpha * 1000, abs=1E-6) == 1.016
    assert drude.thole == 2.6
