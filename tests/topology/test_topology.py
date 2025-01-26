#!/usr/bin/env python3

import os
import pytest
import tempfile
import filecmp
import shutil
from mstk.topology import Topology, UnitCell
from mstk.forcefield import ForceField

cwd = os.path.dirname(os.path.abspath(__file__))


def test_assign_charge():
    top = Topology.open(cwd + '/files/c_3oh.msd')
    ff = ForceField.open(cwd + '/files/c_3oh.ppf')
    ff.assign_charge(top)
    charges = [atom.charge for atom in top.atoms]
    assert pytest.approx(charges, abs=1E-6) == [-0.32, -0.16, 0.4791, -0.45, 0.16, 0.16, 0.16, -0.0291]

    # the function to transfer BCI terms is removed
    # ff.assign_charge(top, transfer_bci_terms=True)
    # charges = [atom.charge for atom in top.atoms]
    # assert pytest.approx(charges, abs=1E-6) == [-0.32, -0.1954, 0.5145, -0.45, 0.16, 0.16, 0.16, -0.0291]


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


def test_packmol():
    # use fixed path so that _pack.inp has fixed content
    tmpdir = os.path.join(tempfile.gettempdir(), '_mstk_' + test_packmol.__name__)
    os.makedirs(tmpdir, exist_ok=True)
    top = Topology.open(cwd + '/files/Im11.zmat')
    top.cell.set_box([3, 3, 3])
    top.scale_with_packmol(10, tempdir=tmpdir)
    assert filecmp.cmp(tmpdir + '/_pack.inp', cwd + '/files/baselines/_pack.inp')
    assert filecmp.cmp(tmpdir + '/_MO_0.xyz', cwd + '/files/baselines/_MO_0.xyz')
    shutil.rmtree(tmpdir)


def test_guess_bonds():
    top = Topology.open(cwd + '/files/MoS2-13x8-layer1.xyz')
    top.cell.set_box([4.109, 4.380, 1.230])
    ff = ForceField.open(cwd + '/files/MoS2.ff')
    top.guess_bonds_from_ff(ff, pbc='xy')
    assert top.n_bond == 1248
