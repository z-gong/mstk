#!/usr/bin/env python3

import os
import tempfile
import filecmp
import sys
import pytest
from mstools.topology import Topology, UnitCell
from mstools.forcefield import ForceField
from mstools.simsys import System

cwd = os.path.dirname(os.path.abspath(__file__))
tmpdir = tempfile.mkdtemp()


def test_gmx():
    ff = ForceField.open(cwd + '/files/10-benzene.ppf')
    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    top.assign_charge_from_ff(ff)
    system = System(top, ff)
    tmpgro = os.path.join(tmpdir, 'conf.gro')
    tmptop = os.path.join(tmpdir, 'topol.top')
    tmpmdp = os.path.join(tmpdir, 'grompp.mdp')
    system.export_gromacs(gro_out=tmpgro, top_out=tmptop, mdp_out=tmpmdp)
    assert filecmp.cmp(tmpgro, cwd + '/files/baselines/conf.gro')
    assert filecmp.cmp(tmptop, cwd + '/files/baselines/topol.top')
    assert filecmp.cmp(tmpmdp, cwd + '/files/baselines/grompp.mdp')


def test_gmx_drude():
    ff = ForceField.open(cwd + '/../forcefield/files/CLP.ff', cwd + '/../forcefield/files/CLPol-alpha.ff')
    top = Topology.open(cwd + '/files/5-Im21-BF4-drude.lmp')
    top.generate_angle_dihedral_improper()
    top.generate_drude_particles(ff)
    top.assign_charge_from_ff(ff)
    system = System(top, ff)
    tmpgro = os.path.join(tmpdir, 'conf-drude.gro')
    tmptop = os.path.join(tmpdir, 'topol-drude.top')
    tmpmdp = os.path.join(tmpdir, 'grompp-drude.mdp')
    system.export_gromacs(gro_out=tmpgro, top_out=tmptop, mdp_out=tmpmdp)
    assert filecmp.cmp(tmpgro, cwd + '/files/baselines/conf-drude.gro')
    assert filecmp.cmp(tmptop, cwd + '/files/baselines/topol-drude.top')
    assert filecmp.cmp(tmpmdp, cwd + '/files/baselines/grompp-drude.mdp')


def test_lmp_drude():
    ff = ForceField.open(cwd + '/../forcefield/files/CLP.ff', cwd + '/../forcefield/files/CLPol-alpha.ff')
    top = Topology.open(cwd + '/files/5-Im21-BF4-drude.lmp')
    top.generate_angle_dihedral_improper()
    top.generate_drude_particles(ff)
    top.assign_charge_from_ff(ff)
    system = System(top, ff)
    tmpdata = os.path.join(tmpdir, 'data-drude.lmp')
    tmpin = os.path.join(tmpdir, 'in-drude.lmp')
    system.export_lammps(data_out=tmpdata, in_out=tmpin)
    assert filecmp.cmp(tmpdata, cwd + '/files/baselines/data-drude.lmp')
    assert filecmp.cmp(tmpin, cwd + '/files/baselines/in-drude.lmp')


def test_charmm():
    ff = ForceField.open(cwd + '/files/10-benzene.ppf')
    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    top.assign_charge_from_ff(ff)
    system = System(top, ff)
    tmppdb = os.path.join(tmpdir, 'conf.pdb')
    tmppsf = os.path.join(tmpdir, 'topol.psf')
    tmpprm = os.path.join(tmpdir, 'ff.prm')
    system.export_charmm(pdb_out=tmppdb, psf_out=tmppsf, prm_out=tmpprm)
    assert filecmp.cmp(tmppdb, cwd + '/files/baselines/conf.pdb')
    assert filecmp.cmp(tmppsf, cwd + '/files/baselines/topol.psf')
    assert filecmp.cmp(tmpprm, cwd + '/files/baselines/ff.prm')
