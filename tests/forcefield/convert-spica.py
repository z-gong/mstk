#!/usr/bin/env python3

import json
from mstk.forcefield import *

ff = ForceField()
ff.vdw_cutoff = 1.5
ff.vdw_long_range = ForceField.VDW_LONGRANGE_SHIFT
ff.lj_mixing_rule = ForceField.LJ_MIXING_NONE
ff.scale_14_vdw = 1.0
ff.scale_14_coulomb = 1.0

with open('spica_top.json') as f:
    j = json.load(f)
for k, v in j['topo'].items():
    for i, typ in enumerate(v['type']):
        if typ in ff.atom_types:
            continue
        atype = AtomType(typ)
        atype.charge = round(v['charge'][i] * math.sqrt(80), 4)
        atype.mass = round(v['mass'][i], 4)
        ff.add_term(atype)
        ### fix charge of COO bead
        if typ == 'COO':
            atype.charge = -1.0
        ##########################

#### add SO4 bead
aterm = AtomType('SO4')
aterm.charge = -1.0
aterm.mass = 96.062
ff.add_term(aterm)
#################

with open('spica_par.json') as f:
    j = json.load(f)
for par in j['params']:
    if par['param'] == 'bond':
        if par['potential'] != 'harmonic':
            raise Exception('Unknown potential form: %s' % par['potential'])
        types = par['types']
        k = par['k'] * 4.184 * 100
        length = par['r0'] / 10
        bterm = HarmonicBondTerm(*types, length, k)
        ff.add_term(bterm)
    elif par['param'] == 'angle':
        types = par['types']
        k = par['k'] * 4.184
        theta = par['theta0'] * DEG2RAD
        if par['potential'] == 'harmonic':
            aterm = HarmonicAngleTerm(*types, theta, k)
        elif par['potential'] == 'sdk':
            aterm = SDKAngleTerm(*types, theta, k)
        else:
            raise Exception('Unknown potential form: %s' % par['potential'])
        ff.add_term(aterm)
    elif par['param'] == 'dihedral':
        if par['potential'] != 'charmm':
            raise Exception('Unknown potential form: %s' % par['potential'])
        types = par['types']
        k = par['k'] * 4.184
        n = par['n']
        phi = par['d'] * DEG2RAD
        dterm = PeriodicDihedralTerm(*types)
        dterm.add_phase(phi, k, n)
        ff.add_term(dterm)
    elif par['param'] == 'pair':
        types = par['types']
        epsilon = par['epsilon'] * 4.184
        sigma = par['sigma'] / 10
        if par['potential'] == 'lj9_6':
            vdw = MieTerm(*types, epsilon, sigma, 9, 6)
        elif par['potential'] == 'lj12_4':
            vdw = MieTerm(*types, epsilon, sigma, 12, 4)
        else:
            raise Exception('Unknown potential form: %s' % par['potential'])
        ff.add_term(vdw)

ff.comments.append('Converted from spica_top.json and spica_prm.json in spica_v1.0.tar.gz')
ff.comments.append('Downloaded from https://www.spica-ff.org/download.html')
ff.write('SPICA_v1.0.zff')
