#!/usr/bin/env python3

import tempfile
import filecmp
import pytest
from mstk.chem.constant import *
from mstk.forcefield import ForceField, Ppf

import os
import shutil

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    ff = ForceField.open(cwd + '/files/TEAM_IL.ppf', cwd + '/files/SPCE.ppf')
    assert len(ff.atom_types) == 65

    c_4pp = ff.atom_types.get('c_4pp')
    assert c_4pp.charge == 0
    assert c_4pp.eqt_vdw == c_4pp.name
    assert c_4pp.eqt_bci == c_4pp.name
    assert c_4pp.eqt_bond == 'c_4'
    assert c_4pp.eqt_ang_c == 'c_4'
    assert c_4pp.eqt_ang_s == 'c_4'
    assert c_4pp.eqt_dih_c == 'c_4'
    assert c_4pp.eqt_dih_s == 'c_4'
    assert c_4pp.eqt_imp_c == 'c_4'
    assert c_4pp.eqt_imp_s == 'c_4'
    assert c_4pp.version == '0.25'

    c_3_ = ff.atom_types.get('c_3-')
    assert c_3_.charge == -1
    assert c_3_.version == '0.1'

    o_1_t = ff.atom_types.get('o_1-t')
    assert o_1_t.charge == -0.3333

    binc = ff.bci_terms.get('c_35an2,h_1')
    assert binc.type1 == 'c_35an2'
    assert binc.type2 == 'h_1'
    assert binc.version == '0.16'
    assert binc.value == -0.22

    vdw = ff.vdw_terms.get('b_4-,b_4-')
    assert vdw.type1 == 'b_4-'
    assert vdw.type2 == 'b_4-'
    assert pytest.approx(vdw.epsilon / 4.184, abs=1E-6) == 0.05
    assert pytest.approx(vdw.sigma * 2 ** (1 / 6) * 10, abs=1E-6) == 3.8
    assert vdw.version == '0.36'

    vdw = ff.pairwise_vdw_terms.get('o_2,o_2w')
    assert vdw.type1 == 'o_2'
    assert vdw.type2 == 'o_2w'
    assert pytest.approx(vdw.epsilon / 4.184, abs=1E-6) == 0.22936
    assert pytest.approx(vdw.sigma * 2 ** (1 / 6) * 10, abs=1E-6) == 3.44671
    assert vdw.version == None

    bond = ff.bond_terms.get('c_2n,n_2-')
    assert bond.type1 == 'c_2n'
    assert bond.type2 == 'n_2-'
    assert pytest.approx(bond.length * 10, abs=1E-6) == 1.276
    assert pytest.approx(bond.k / 4.184 / 100, abs=1E-6) == 487.36
    assert bond.version == '0.3'

    angle = ff.angle_terms.get('c_3a,c_3a,h_1')
    assert angle.type1 == 'c_3a'
    assert angle.type2 == 'c_3a'
    assert angle.type3 == 'h_1'
    assert pytest.approx(angle.theta, abs=1E-6) == 118.916 * DEG2RAD
    assert pytest.approx(angle.k / 4.184, abs=1E-6) == 43.8769
    assert angle.version == '0.1'

    dihedral = ff.dihedral_terms.get('*,c_4,c_4,*')
    assert dihedral.type1 == '*'
    assert dihedral.type2 == 'c_4'
    assert dihedral.type3 == 'c_4'
    assert dihedral.type4 == '*'
    phase1 = dihedral.phases[0]
    assert pytest.approx(phase1.k / 4.184, abs=1E-6) == 0.6214
    assert phase1.n == 1
    assert phase1.phi == 0.0
    phase2 = dihedral.phases[1]
    assert pytest.approx(phase2.k / 4.184, abs=1E-6) == 0.0507
    assert phase2.n == 2
    assert pytest.approx(phase2.phi, abs=1E-6) == PI
    phase3 = dihedral.phases[2]
    assert pytest.approx(phase3.k / 4.184, abs=1E-6) == 0.0688
    assert phase3.n == 3
    assert phase3.phi == 0.0
    assert dihedral.version == '0.12'

    improper = ff.improper_terms.get('c_3o,o_1,*,*')
    assert pytest.approx(improper.k / 4.184, abs=1E-6) == 18
    assert improper.version == '0.21'

    assert len(ff.polarizable_terms) == 0
