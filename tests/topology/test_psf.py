#!/usr/bin/env python3

import os
import tempfile
import filecmp
import pytest
from mstools.topology import Topology

cwd = os.path.dirname(os.path.abspath(__file__))
tmpdir = tempfile.mkdtemp()


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


def test_write():
    lmp = Topology.open(cwd + '/files/10-H2O-5-C3H6.lmp', improper_center=3)
    tmp = os.path.join(tmpdir, 'lmp-out.psf')
    lmp.write(tmp)
    assert filecmp.cmp(tmp, cwd + '/files/baselines/lmp-out.psf')

    zmat = Topology.open(cwd + '/files/Im11.zmat')
    tmp = os.path.join(tmpdir, 'zmat-out.psf')
    zmat.write(tmp)
    assert filecmp.cmp(tmp, cwd + '/files/baselines/zmat-out.psf')
