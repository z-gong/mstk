#!/usr/bin/env python3

import pytest
from mstk.chem.constant import *
from mstk.forcefield import ForceField, Padua, PaduaLJScaler

import os

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    ff = ForceField.open(cwd + '/files/CLP.ff', cwd + '/files/CLPol-alpha.ff')
    br = ff.atom_types['Br']
    assert br.mass == 79.904
    assert br.charge == -1.0

    cap = ff.atom_types['CAP']
    ca = ff.atom_types['CA']
    assert cap.eqt_vdw == 'CAP'
    assert cap.eqt_bond == 'CA'
    assert cap.eqt_imp_s == ca.name
    assert ca.eqt_vdw == ca.name
    assert ca.eqt_ang_c == ca.name
    assert ca.eqt_dih_c == 'CA'

    vdw = ff.vdw_terms['Br,Br']
    assert vdw.epsilon == 0.86
    assert vdw.sigma == 3.97 / 10

    bond = ff.bond_terms['CT,CT']
    assert bond.length == 1.529 / 10
    assert pytest.approx(bond.k, abs=1E-6) == 2242.0 / 2 * 100

    angle = ff.angle_terms['CT,CT,HC']
    assert pytest.approx(angle.theta, abs=1e-6) == 110.7 * DEG2RAD
    assert angle.k == 313.8 / 2

    term = ff.dihedral_terms['CT,CT,CT,HC']
    assert term.k1 == 0
    assert term.k2 == 0
    assert term.k3 == 1.2552 / 2
    assert term.k4 == 0

    term = ff.dihedral_terms['CT,CT,CT,CT']
    assert term.k1 == 5.4392 / 2
    assert term.k2 == -0.2092 / 2
    assert term.k3 == 0.8368 / 2
    assert term.k4 == 0

    improper = ff.improper_terms['NA,CR,CT,CW']
    assert improper.k == 8.368 / 2

    drude = ff.polar_terms['CT']
    assert drude.mass == 0.4
    assert pytest.approx(drude.k / 100, abs=1E-6) == 4184 / 2
    assert pytest.approx(drude.alpha * 1000, abs=1E-6) == 1.016
    assert drude.thole == 2.6


def test_ljscaler():
    ff = ForceField.open(cwd + '/files/CLP.ff', cwd + '/files/CLPol-alpha.ff')
    scaler = PaduaLJScaler(cwd + '/files/CLPol-ljscale.ff')
    scaler.scale(ff)

    vdw = ff.get_vdw_term('C1', 'C1')
    assert pytest.approx(vdw.epsilon, abs=1E-5) == 0.047123 * 4.184
    assert pytest.approx(vdw.sigma, abs=1E-5) == 3.869688 / 10 / 2 ** (1 / 6)
    assert vdw.comments == ['eps*0.714', 'sig*0.985']

    vdw = ff.get_vdw_term('CR', 'CZA')
    assert pytest.approx(vdw.epsilon, abs=1E-5) == 0.046627 * 4.184
    assert pytest.approx(vdw.sigma, abs=1E-5) == 3.784243 / 10 / 2 ** (1 / 6)
    assert vdw.comments == ['eps*0.686', 'sig*0.985']

