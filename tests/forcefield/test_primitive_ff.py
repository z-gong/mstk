#!/usr/bin/env python3

import os
from mstools.topology import Topology, Molecule, Atom, Bond
from mstools.forcefield import ForceField
from mstools.forcefield.typer import ZftTyper, typer_primitive
from mstools.simsys import System

cwd = os.path.dirname(os.path.abspath(__file__))


def test_primitive_ff():
    mol = Molecule.from_smiles('C=CCC#CCOC(=O)NC(O)C(=O)OC')
    typer_primitive.type(mol)
    assert [atom.type for atom in mol.atoms] == [
        'C3', 'C3', 'C4', 'C2', 'C2', 'C4', 'O2', 'C3=o', 'O1', 'N3co', 'C4', 'O2', 'C3=o', 'O1', 'O2', 'C4',
        'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1po', 'H1', 'H1po', 'H1', 'H1', 'H1'
    ]
    top = Topology([mol])
    ff = ForceField.open('primitive.zff')
    ff.assign_charge(top)
    system = System(top, ff)
