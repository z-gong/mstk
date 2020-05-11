#!/usr/bin/env python3

from mstools.forcefield import Padua, Ppf, Zfp, ZftTyper
from mstools.forcefield.ffterm import *
from mstools.topology import *

import os

cwd = os.path.dirname(os.path.abspath(__file__))


def test_typing():
    im61 = Molecule.from_smiles('C[n+]1cn(cc1)CCCCCC')
    im2eben = Molecule.from_smiles('C[n+]1cn(cc1)CCc1ccccc1')
    im2bben = Molecule.from_smiles('C[n+]1cn(cc1)CCCc1ccccc1')
    dca = Molecule.from_smiles('N#C[N-]C#N')
    fsi = Molecule.from_smiles('FS(=O)(=O)[N-]S(=O)(=O)F')
    tfsi = Molecule.from_smiles('FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)(F)')
    top = Topology([im61, im2eben, im2bben, dca, fsi, tfsi])

    typer = ZftTyper(cwd + '/files/CLP-define.zft')
    typer.type(top)

    top.write(cwd + '/files/zft.xyz')

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
