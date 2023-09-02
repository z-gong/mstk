#!/usr/bin/env python3

import os
import tempfile
import filecmp
import pytest
from mstk.topology import Topology, UnitCell
from mstk.forcefield import ForceField
from mstk.simsys import System
import shutil

cwd = os.path.dirname(os.path.abspath(__file__))


def test_gmx():
    tmpdir = tempfile.mkdtemp()
    ff = ForceField.open(cwd + '/files/10-benzene.ppf')
    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    ff.assign_charge(top)
    system = System(top, ff)
    tmpgro = os.path.join(tmpdir, 'conf.gro')
    tmptop = os.path.join(tmpdir, 'topol.top')
    tmpmdp = os.path.join(tmpdir, 'grompp.mdp')
    system.export_gromacs(gro_out=tmpgro, top_out=tmptop, mdp_out=tmpmdp)
    assert filecmp.cmp(tmpgro, cwd + '/files/baselines/conf.gro')
    assert filecmp.cmp(tmptop, cwd + '/files/baselines/topol.top')
    assert filecmp.cmp(tmpmdp, cwd + '/files/baselines/grompp.mdp')
    shutil.rmtree(tmpdir)


def test_gmx_drude():
    tmpdir = tempfile.mkdtemp()
    ff = ForceField.open(cwd + '/../forcefield/files/CLP.ff', cwd + '/../forcefield/files/CLPol-alpha.ff')
    top = Topology.open(cwd + '/files/5-Im21-BF4-drude.lmp')
    top.generate_angle_dihedral_improper()
    top.generate_drude_particles(ff)
    ff.assign_charge(top)
    system = System(top, ff)
    tmpgro = os.path.join(tmpdir, 'conf-drude.gro')
    tmptop = os.path.join(tmpdir, 'topol-drude.top')
    tmpmdp = os.path.join(tmpdir, 'grompp-drude.mdp')
    system.export_gromacs(gro_out=tmpgro, top_out=tmptop, mdp_out=tmpmdp)
    assert filecmp.cmp(tmpgro, cwd + '/files/baselines/conf-drude.gro')
    assert filecmp.cmp(tmptop, cwd + '/files/baselines/topol-drude.top')
    assert filecmp.cmp(tmpmdp, cwd + '/files/baselines/grompp-drude.mdp')
    shutil.rmtree(tmpdir)


def test_gmx_tip4p():
    tmpdir = tempfile.mkdtemp()
    mol = Topology.open(cwd + '/../topology/files/TIP3P.zmat').molecules[0]
    ff = ForceField.open(cwd + '/../forcefield/files/TIP4P.zfp')
    mol.generate_virtual_sites(ff)
    ff.assign_charge(mol)

    top = Topology([mol], numbers=[100], cell=UnitCell([3, 3, 3]))
    top.set_positions(Topology.open(cwd + '/files/100-TIP4P.pdb').positions)

    system = System(top, ff)
    tmpgro = os.path.join(tmpdir, 'conf-vsite.gro')
    tmptop = os.path.join(tmpdir, 'topol-vsite.top')
    tmpmdp = os.path.join(tmpdir, 'grompp-vsite.mdp')
    system.export_gromacs(gro_out=tmpgro, top_out=tmptop, mdp_out=tmpmdp)
    assert filecmp.cmp(tmpgro, cwd + '/files/baselines/conf-vsite.gro')
    assert filecmp.cmp(tmptop, cwd + '/files/baselines/topol-vsite.top')
    assert filecmp.cmp(tmpmdp, cwd + '/files/baselines/grompp-vsite.mdp')
    shutil.rmtree(tmpdir)


def test_namd():
    tmpdir = tempfile.mkdtemp()
    ff = ForceField.open(cwd + '/files/10-benzene.ppf')
    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    ff.assign_charge(top)
    system = System(top, ff)
    tmppdb = os.path.join(tmpdir, 'conf.pdb')
    tmppsf = os.path.join(tmpdir, 'topol.psf')
    tmpprm = os.path.join(tmpdir, 'ff.prm')
    system.export_namd(pdb_out=tmppdb, psf_out=tmppsf, prm_out=tmpprm)
    assert filecmp.cmp(tmppdb, cwd + '/files/baselines/conf.pdb')
    assert filecmp.cmp(tmppsf, cwd + '/files/baselines/topol.psf')
    assert filecmp.cmp(tmpprm, cwd + '/files/baselines/ff.prm')
    shutil.rmtree(tmpdir)
