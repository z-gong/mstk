#!/usr/bin/env python3

import os
import pytest
import tempfile
import filecmp
import shutil
from mstools.topology import Topology, UnitCell
from mstools.forcefield import ForceField

cwd = os.path.dirname(os.path.abspath(__file__))


def test_assign_charge():
    top = Topology.open(cwd + '/files/c_3oh.msd')
    ff = ForceField.open(cwd + '/files/c_3oh.ppf')
    top.assign_charge_from_ff(ff)
    charges = [atom.charge for atom in top.atoms]
    assert pytest.approx(charges, abs=1E-6) == [-0.32, -0.16, 0.4791, -0.45, 0.16, 0.16, 0.16, -0.0291]

    top.assign_charge_from_ff(ff, transfer_qinc_terms=True)
    charges = [atom.charge for atom in top.atoms]
    assert pytest.approx(charges, abs=1E-6) == [-0.32, -0.1954, 0.5145, -0.45, 0.16, 0.16, 0.16, -0.0291]


def test_compress():
    top = Topology.open(cwd + '/files/10-H2O-5-C3H6.lmp', improper_center=3)
    molecules = top.molecules
    for mol in molecules[4:]:
        for atom in mol.atoms:
            atom.charge *= 2

    top = Topology.open(cwd + '/files/Im11.zmat')
    molecules += top.molecules

    top = Topology.open(cwd + '/files/10-H2O-5-C3H6.lmp', improper_center=3)
    molecules += top.molecules

    top = Topology.open(cwd + '/files/Im11.zmat')
    molecules += top.molecules

    top = Topology()
    top.update_molecules(molecules)
    mols_unique = top.get_unique_molecules()
    for mol, count in mols_unique.items():
        print(str(mol), count)

    assert list(mols_unique.values()) == [4, 6, 5, 1, 10, 5, 1]


def test_scale():
    tmpdir = tempfile.mkdtemp()
    top = Topology.open(cwd + '/files/Im11.zmat')
    top.cell.set_box([3, 3, 3])
    os.chdir(tmpdir)
    top.scale_with_packmol(10)
    os.chdir(cwd)
    assert filecmp.cmp(tmpdir + '/_pack.inp', cwd + '/files/baselines/_pack.inp')
    assert filecmp.cmp(tmpdir + '/_MO_0.xyz', cwd + '/files/baselines/_MO_0.xyz')
    shutil.rmtree(tmpdir)


def test_guess_connectivity():
    top = Topology.open(cwd + '/files/MoS2-13x8-layer1.xyz')
    top.cell.set_box([4.109, 4.380, 1.230])
    ff = ForceField.open(cwd + '/files/MoS2.ff')
    top.guess_connectivity_from_ff(ff, angle_tolerance=15, pbc='xy')
    assert top.n_bond == 1248
    assert top.n_angle == 3120
    assert top.n_dihedral == 0
    assert top.n_improper == 0
