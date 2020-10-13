#!/usr/bin/env python3

import os
import tempfile
import filecmp
import pytest
from mstools.topology import Topology, Molecule

cwd = os.path.dirname(os.path.abspath(__file__))


def test_smiles():
    im41 = Molecule.from_smiles('C[n+]1cn(cc1)CCCC Im41')
    assert im41.name == 'Im41'
    assert im41.n_atom == 25
    assert im41.n_improper == 5

    bf4 = Molecule.from_smiles('[B-](F)(F)(F)F')
    assert bf4.n_atom == 5
    assert bf4.name == 'BF4-'
    assert bf4.n_improper == 0
