#!/usr/bin/env python3

import math
from mstk.forcefield import ForceField
from mstk.forcefield.ffterm import *


def main():
    with open('gaff.dat') as f:
        lines = f.read().splitlines()
    lines_atom = lines[1:84]
    lines_bond = lines[86:1014]
    lines_angle = lines[1015:6330]
    lines_dihedral = lines[6331:7075]
    lines_improper = lines[7076:7114]
    lines_vdw = lines[7119:7202]

    ff = ForceField()
    ff.vdw_cutoff = 1.0
    ff.vdw_long_range = ff.VDW_LONGRANGE_CORRECT
    ff.lj_mixing_rule = ff.LJ_MIXING_LB
    ff.scale_14_vdw = 0.5
    ff.scale_14_coulomb = 0.8333

    for line in lines_atom:
        atype, str_mass = line.split()[:2]
        aterm = AtomType(atype, mass=float(str_mass))
        ff.add_term(aterm)
    for line in lines_bond:
        atypes = line[:2].strip(), line[3:5].strip()
        str_k, str_length = line[5:].split()[:2]
        k = float(str_k) * 100 * 4.184
        length = float(str_length) / 10
        bterm = HarmonicBondTerm(*atypes, length, k)
        ff.add_term(bterm)
    for line in lines_angle:
        atypes = line[:2].strip(), line[3:5].strip(), line[6:8].strip()
        str_k, str_theta = line[8:].split()[:2]
        k = float(str_k) * 4.184
        theta = float(str_theta) * DEG2RAD
        aterm = HarmonicAngleTerm(*atypes, theta, k)
        ff.add_term(aterm)
    for line in lines_dihedral:
        type1, type2, type3, type4 = line[:2].strip(), line[3:5].strip(), line[6:8].strip(), line[9:11].strip()
        str_divider, str_k, str_phi, str_n = line[11:].split()[:4]
        k = float(str_k) / int(str_divider) * 4.184
        phi = float(str_phi) * DEG2RAD
        n = abs(int(float(str_n)))
        dterm = PeriodicDihedralTerm(type1.replace('X', '*'), type2, type3, type4.replace('X', '*'))
        if dterm.name in ff.dihedral_terms and k == 0:
            continue
        dterm = ff.dihedral_terms.get(dterm.name, dterm)
        dterm.add_parameter(phi, k, n)
        ff.add_term(dterm, replace=True)
    for line in lines_improper:
        type1, type2, type3, type4 = line[:2].strip(), line[3:5].strip(), line[6:8].strip(), line[9:11].strip()
        str_k = line[11:].split()[0]
        k = float(str_k) * 4.184
        iterm = OplsImproperTerm(type3, type1.replace('X', '*'), type2.replace('X', '*'), type4.replace('X', '*'), k)
        ff.add_term(iterm, replace=True)
    for line in lines_vdw:
        atype, str_r0_2, str_epsilon = line.split()[:3]
        sigma = float(str_r0_2) * 2 / 2 ** (1 / 6) / 10
        epsilon = float(str_epsilon) * 4.184
        ljterm = LJ126Term(atype, atype, epsilon, sigma)
        ff.add_term(ljterm)

    ff.write('gaff.zfp')


if __name__ == '__main__':
    main()
