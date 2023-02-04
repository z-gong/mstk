#!/usr/bin/env python3

import os
from mstk.topology import Topology, Molecule, Atom, Bond
from mstk.forcefield import ForceField
from mstk.forcefield.typer import ZftTyper, typer_primitive
from mstk.simsys import System

cwd = os.path.dirname(os.path.abspath(__file__))


def test_primitive_ff():
    mol = Molecule.from_smiles('C=CCC#CCOC(=O)NC(O)C(=O)OC')
    typer_primitive.type(mol)
    assert [atom.type for atom in mol.atoms] == [
        'C3', 'C3', 'C4', 'C2', 'C2', 'C4', 'O2', 'C3=O', 'O1', 'N3CO', 'C4', 'O2', 'C3=O', 'O1', 'O2', 'C4',
        'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1p', 'H1', 'H1p', 'H1', 'H1', 'H1'
    ]
    top = Topology([mol])
    ff = ForceField.open('primitive.zff')
    ff.assign_charge(top)
    system = System(top, ff)
