#!/usr/bin/env python3

import os
import tempfile
import filecmp
import pytest
import shutil
from mstk.topology import Topology, TIP4PSite
from mstk.forcefield import ForceField

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    psf = Topology.open(cwd + '/files/10-H2O-5-C3H6.psf')
    assert psf.is_drude == False
    assert psf.n_atom == 75
    assert psf.n_molecule == 15
    assert psf.n_bond == 60
    assert psf.n_angle == 70
    assert psf.n_dihedral == 50
    assert psf.n_improper == 10
    assert psf.cell.volume == 0

    assert psf.bonds[-1].atom1.name == 'C69'
    assert psf.angles[-1].atom2.name == 'C69'
    assert psf.dihedrals[-1].atom4.name == 'H75'
    assert psf.impropers[-1].atom1.type == 'c_3h'

    atom = psf.atoms[0]
    assert atom.id == 0
    assert atom.molecule.id == 0
    assert atom.type == 'o_2w'
    assert atom.name == 'O1'
    assert atom.charge == -0.8476
    assert atom.mass == 15.9994

    atom = psf.atoms[-1]
    assert atom.id == 74
    assert atom.molecule.id == 14
    assert atom.type == 'h_1'
    assert atom.name == 'H75'
    assert atom.charge == 0.06
    assert atom.mass == 1.00794

    mol = psf.molecules[0]
    assert mol.id == 0
    assert mol.name == 'WAT'

    mol = psf.molecules[-1]
    assert mol.id == 14
    assert mol.name == 'C3H'


def test_read_swm4():
    psf = Topology.open(cwd + '/files/9-SWM4.psf')
    assert psf.n_molecule == 9

    mol = psf.molecules[-1]
    assert mol.n_atom == 5

    assert not mol.atoms[0].is_drude
    assert mol.atoms[1].is_drude
    assert not mol.atoms[2].is_drude
    assert not mol.atoms[3].is_drude
    assert not mol.atoms[4].is_drude

    atom = mol.atoms[0]
    assert atom.alpha == 0.0009782
    assert atom.thole == 2.6

    assert mol.atoms[0].virtual_site is None
    assert mol.atoms[1].virtual_site is None
    assert mol.atoms[2].virtual_site is None
    assert mol.atoms[3].virtual_site is None
    assert mol.atoms[4].virtual_site is not None

    vsite = mol.atoms[-1].virtual_site
    assert type(vsite) is TIP4PSite
    assert vsite.parents[0] is mol.atoms[0]
    assert vsite.parents[1] is mol.atoms[2]
    assert vsite.parents[2] is mol.atoms[3]
    assert vsite.parameters == [0.024034]


def test_write():
    tmpdir = tempfile.mkdtemp()
    lmp = Topology.open(cwd + '/files/10-H2O-5-C3H6.lmp', improper_center=3)
    tmp = os.path.join(tmpdir, 'lmp-out.psf')
    lmp.write(tmp)
    assert filecmp.cmp(tmp, cwd + '/files/baselines/lmp-out.psf')

    zmat = Topology.open(cwd + '/files/Im11.zmat')
    tmp = os.path.join(tmpdir, 'zmat-out.psf')
    zmat.write(tmp)
    assert filecmp.cmp(tmp, cwd + '/files/baselines/zmat-out.psf')
    shutil.rmtree(tmpdir)


def test_write_swm4():
    tmpdir = tempfile.mkdtemp()
    ff = ForceField.open(cwd + '/../forcefield/files/SWM4-NDP.zfp')
    mol = Topology.open(cwd + '/files/TIP3P.zmat').molecules[0]
    mol.generate_virtual_sites(ff)
    mol.generate_drude_particles(ff)
    mol.assign_charge_from_ff(ff)

    top = Topology([mol], numbers=[10])
    tmp = os.path.join(tmpdir, 'swm4-out.psf')
    top.write(tmp)

    assert filecmp.cmp(tmp, cwd + '/files/baselines/swm4-out.psf')
    shutil.rmtree(tmpdir)
