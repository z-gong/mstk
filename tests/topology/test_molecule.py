#!/usr/bin/env python3

import os
import math
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


def test_connectivity():
    ethane = Topology.open(cwd + '/files/CH3NH2.pdb').molecules[0]
    assert pytest.approx(ethane.bonds[0].evaluate(), abs=1E-4) == 0.1070
    assert pytest.approx(ethane.angles[0].evaluate(), abs=1E-4) == 109.5278 / 180 * math.pi
    assert pytest.approx(ethane.dihedrals[0].evaluate(), abs=1E-4) == 60.0301 / 180 * math.pi
    assert pytest.approx(ethane.dihedrals[-1].evaluate(), abs=1E-4) == -60.0159 / 180 * math.pi
    assert pytest.approx(ethane.impropers[0].evaluate(), abs=1E-4) == -33.0905 / 180 * math.pi


def test_distance_matrix():
    ethanol = Molecule.from_smiles('CO ethanol')
    assert ethanol.get_distance_matrix().tolist() == [[0, 1, 1, 1, 1, 2],
                                                      [1, 0, 2, 2, 2, 1],
                                                      [1, 2, 0, 2, 2, 3],
                                                      [1, 2, 2, 0, 2, 3],
                                                      [1, 2, 2, 2, 0, 3],
                                                      [2, 1, 3, 3, 3, 0]]
    assert ethanol.get_distance_matrix(max_bond=2).tolist() == [[0, 1, 1, 1, 1, 2],
                                                                [1, 0, 2, 2, 2, 1],
                                                                [1, 2, 0, 2, 2, 0],
                                                                [1, 2, 2, 0, 2, 0],
                                                                [1, 2, 2, 2, 0, 0],
                                                                [2, 1, 0, 0, 0, 0]]
