#!/usr/bin/env python3

import os
import sys
import pytest
from mstools.topology import Topology, UnitCell
from mstools.forcefield import Ppf, Padua, Zfp
from mstools.simsys import System

import simtk.openmm as mm
from simtk.openmm import app
from simtk.unit import kilocalorie_per_mole as kcal_mol

cwd = os.path.dirname(os.path.abspath(__file__))


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
    top.init_from_molecules(molecules)
    mols_unique = top.get_unique_molecules()
    for mol, count in mols_unique.items():
        print(str(mol), count)

    assert list(mols_unique.values()) == [4, 6, 5, 1, 10, 5, 1]
