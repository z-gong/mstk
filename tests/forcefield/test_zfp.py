#!/usr/bin/env python3

import os
import pytest

from mstools.forcefield import ForceField, Zfp
from mstools.forcefield.ffterm import *

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    ff = ForceField.open(cwd + '/files/TEAM_IL.zfp')
    assert ff.vdw_cutoff == 1.2
    assert ff.vdw_long_range == 'correct'
    assert ff.scale_14_vdw == 0.5
    assert ff.lj_mixing_rule == ForceField.LJ_MIXING_LB
    assert pytest.approx(ff.scale_14_coulomb, abs=1E-4) == 0.8333

    atype = ff.atom_types['s_1-']
    assert atype.name == 's_1-'
    assert atype.charge == -1
    assert atype.eqt_imp_s == 's_1'
    assert atype.eqt_imp_c == 's_1-'

    term = ff.charge_increment_terms['b_4-,c_2nb']
    assert term.type1 == 'b_4-'
    assert term.type2 == 'c_2nb'
    assert term.value == 0.262

    term = ff.charge_increment_terms['b_4-,f_1']
    assert term.type1 == 'b_4-'
    assert term.type2 == 'f_1'
    assert term.value == 0.4215

    term = ff.dihedral_terms['o_1,s_4o,n_2-,s_4']
    assert type(term) == PeriodicDihedralTerm
    assert term.type1 == 'o_1'
    assert term.type2 == 's_4o'
    assert term.type3 == 'n_2-'
    assert term.type4 == 's_4'
    assert term.follow_opls_convention
    assert pytest.approx(term.get_opls_parameters(), abs=1E-4) == [12.9809, 7.9709, 1.6389, 0]

    term = ff.improper_terms['c_3a,c_3a,c_3a,h_1']
    assert type(term) == OplsImproperTerm
    assert term.type1 == 'c_3a'
    assert term.type2 == 'c_3a'
    assert term.type3 == 'c_3a'
    assert term.type4 == 'h_1'
    assert pytest.approx(term.k, abs=1E-4) == 7.1270

    ff.write(cwd + '/files/zfp-out.zfp')


def test_write():
    clp = ForceField.open(cwd + '/files/CLP.ff', cwd + '/files/CLPol-alpha.ff')
    Zfp.save_to(clp, cwd + '/files/out-CLPol.zfp')

    il = ForceField.open(cwd + '/files/TEAM_IL.ppf')
    Zfp.save_to(il, cwd + '/files/out-TEAM_IL.zfp')

    spce = ForceField.open(cwd + '/files/SPCE.ppf')
    spce.write(cwd + '/files/out-SPCE.zfp')
