#!/usr/bin/env python3

import os
import math
import pytest
from mstools.topology import Topology, Molecule

cwd = os.path.dirname(os.path.abspath(__file__))


def test_cccc():
    cccc = Molecule.from_smiles('C=CC=C')
    assert [b.is_aromatic for b in cccc.bonds] == [False] * 9
    assert [b.is_conjugated for b in cccc.bonds] == [True] * 3 + [False] * 6


def test_cccco():
    cccc = Molecule.from_smiles('C=CC=CO')
    assert [b.is_aromatic for b in cccc.bonds] == [False] * 10
    assert [b.is_conjugated for b in cccc.bonds] == [True] * 4 + [False] * 6


def test_im1():
    im1 = Molecule.from_smiles('C[n+]1c[nH](cc1)')
    assert [b.is_aromatic for b in im1.bonds] == [False] + [True] * 5 + [False] * 7
    assert [b.is_conjugated for b in im1.bonds] == [False] + [True] * 5 + [False] * 7


def test_bf4():
    bf4 = Molecule.from_smiles('[B-](F)(F)(F)F')
    assert [b.is_aromatic for b in bf4.bonds] == [False] * 4
    assert [b.is_conjugated for b in bf4.bonds] == [False] * 4
