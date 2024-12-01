#!/usr/bin/env python3

from mstk.forcefield.typer import Typer, SmartsTyper
from mstk.topology import Topology, Molecule, Atom, Bond
import os

cwd = os.path.dirname(os.path.abspath(__file__))

typer_primitive = Typer.open('primitive.smt')

def test_smarts_typer():
    im61 = Molecule.from_smiles('C[n+]1cn(cc1)CCCCCC')
    im2eben = Molecule.from_smiles('C[n+]1cn(cc1)CCc1ccccc1')
    im2bben = Molecule.from_smiles('C[n+]1cn(cc1)CCCc1ccccc1')
    dca = Molecule.from_smiles('N#C[N-]C#N')
    fsi = Molecule.from_smiles('FS(=O)(=O)[N-]S(=O)(=O)F')
    tfsi = Molecule.from_smiles('FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)(F)')
    top = Topology([im61, im2eben, im2bben, dca, fsi, tfsi])

    typer = SmartsTyper(cwd + '/files/CLP-define.smt')
    typer.type(top)

    assert [atom.type for atom in im61.atoms] == (
            ['C1', 'NA', 'CR', 'NA', 'CW', 'CW', 'C1', 'C2', 'CS', 'CS', 'CS', 'CT']
            + ['H1', 'H1', 'H1', 'HCR', 'HCW', 'HCW'] + ['H1', 'H1'] + ['HC'] * 11)
    assert [atom.type for atom in im2eben.atoms] == (
            ['C1', 'NA', 'CR', 'NA', 'CW', 'CW', 'C_4H2_NA',
             'C_4H2_CA', 'CAT', 'CAO', 'CAM', 'CAP', 'CAM', 'CAO']
            + ['H1', 'H1', 'H1', 'HCR', 'HCW', 'HCW']
            + ['H1', 'H1', 'HT', 'HT'] + ['HAT'] * 5)
    assert [atom.type for atom in dca.atoms] == ['NZA', 'CZA', 'N3A', 'CZA', 'NZA']
    assert [atom.type for atom in fsi.atoms] == (
        ['FSI', 'SBT', 'OBT', 'OBT', 'NBT', 'SBT', 'OBT', 'OBT', 'FSI'])
    assert [atom.type for atom in tfsi.atoms] == (
        ['F1', 'CBT', 'F1', 'F1', 'SBT', 'OBT', 'OBT', 'NBT', 'SBT', 'OBT', 'OBT',
         'CBT', 'F1', 'F1', 'F1'])


def test_typer_primitive():
    mol = Molecule()
    for i in range(6):
        atom = Atom()
        atom.symbol = 'C'
        mol.add_atom(atom)
    for i in range(6):
        C = mol.atoms[i]
        C2 = mol.atoms[i + 1] if i < 5 else mol.atoms[0]
        mol.add_bond(C, C2, order=Bond.Order.DOUBLE)

    typer_primitive.type(mol)
    assert [atom.type for atom in mol.atoms] == ['C2'] * 6


def test_typer_primitive_kekulize():
    mol = Molecule()
    for i in range(6):
        atom = Atom()
        atom.symbol = 'C'
        mol.add_atom(atom)
    for i in range(6):
        atom = Atom()
        atom.symbol = 'H'
        mol.add_atom(atom)
    for i in range(6):
        C, H = mol.atoms[i], mol.atoms[i + 6]
        mol.add_bond(C, H, order=Bond.Order.SINGLE)
        C2 = mol.atoms[i + 1] if i < 5 else mol.atoms[0]
        order = Bond.Order.SINGLE if i % 2 else Bond.Order.DOUBLE
        mol.add_bond(C, C2, order=order)

    typer_primitive.type(mol)
    assert [atom.type for atom in mol.atoms] == ['C3a'] * 6 + ['H1'] * 6


