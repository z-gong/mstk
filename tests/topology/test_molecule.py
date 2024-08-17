#!/usr/bin/env python3

import os
import pytest
from mstk.chem.constant import *
from mstk.topology import Topology, Molecule, Bond

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
    assert pytest.approx(ethane.angles[0].evaluate(), abs=1E-4) == 109.5278 * DEG2RAD
    assert pytest.approx(ethane.dihedrals[0].evaluate(), abs=1E-4) == 60.0301 * DEG2RAD
    assert pytest.approx(ethane.dihedrals[-1].evaluate(), abs=1E-4) == -60.0159 * DEG2RAD
    assert pytest.approx(ethane.impropers[0].evaluate(), abs=1E-4) == -33.0905 * DEG2RAD


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


def test_merge_split():
    def print_mol(mol):
        print(mol)
        for atom in mol.atoms:
            print(atom.residue, atom)
        for bond in mol.bonds:
            print(bond)

    ethanol = Molecule.from_smiles('CO ethanol')
    water = Molecule.from_smiles('O water')
    ethyne = Molecule.from_smiles('C#C ethyne')
    mol = Molecule.merge([ethanol, water, ethyne])
    assert mol.n_atom == 13

    C1_ethyne = mol.atoms[-4]
    C2_ethyne = mol.atoms[-3]
    assert C2_ethyne.name == 'C2'
    bond = C2_ethyne.bonds[0]
    assert bond.equals(Bond(C1_ethyne, C2_ethyne))
    mol.remove_connectivity(bond)
    mol.generate_angle_dihedral_improper()

    O_ethanol = mol.atoms[1]
    assert O_ethanol.name == 'O2'
    mol.add_bond(O_ethanol, C1_ethyne)

    pieces = mol.split()
    assert len(pieces) == 3
    piece1, piece2, piece3 = pieces
    assert piece1.n_atom == 8
    assert piece1.n_residue == 2
    assert piece1.n_bond == 7
    assert piece1.n_angle == 7
    print_mol(piece1)
    assert piece2.n_atom == 3
    assert piece2.n_residue == 1
    assert piece2.n_bond == 2
    assert piece2.n_angle == 1
    print_mol(piece2)
    assert piece3.n_atom == 2
    assert piece3.n_residue == 1
    assert piece3.n_bond == 1
    assert piece3.n_angle == 0
    print_mol(piece3)

    pieces = mol.split_residues()
    assert len(pieces) == 2
    piece1, piece2 = pieces
    assert piece1.n_atom == 10
    assert piece1.n_residue == 2
    assert piece1.n_bond == 8
    assert piece1.n_angle == 7
    print_mol(piece1)
    assert piece2.n_atom == 3
    assert piece2.n_residue == 1
    assert piece2.n_bond == 2
    assert piece2.n_angle == 1
    print_mol(piece2)