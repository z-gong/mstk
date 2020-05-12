#!/usr/bin/env python3

import os
import sys
import pytest
from mstools.topology import Topology, UnitCell
from mstools.forcefield import Ppf, Padua, Zfp
from mstools.simsys import System

cwd = os.path.dirname(os.path.abspath(__file__))


def test_gmx():
    ff = Ppf(cwd + '/files/10-benzene.ppf')
    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    top.assign_charge_from_ff(ff)
    system = System(top, ff)
    system.export_gromacs(gro_out=cwd + '/files/conf.gro',
                          top_out=cwd + '/files/topol.top',
                          mdp_out=cwd + '/files/grompp.mdp')


def test_gmx_drude():
    pass


def test_charmm():
    ff = Ppf(cwd + '/files/10-benzene.ppf')
    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    top.assign_charge_from_ff(ff)
    system = System(top, ff)
    system.export_charmm(pdb_out=cwd + '/files/conf.pdb',
                         psf_out=cwd + '/files/psf.psf',
                         prm_out=cwd + '/files/ff.prm')
