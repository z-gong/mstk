import os
import tempfile
import filecmp
import pytest

from mstk.topology import Molecule, Topology
from mstk.forcefield import ForceField, Zfp
from mstk.forcefield.typer import ZftTyper
from mstk.simsys import System

cwd = os.path.dirname(os.path.abspath(__file__))


@pytest.mark.parametrize('name_smiles', [
    ('1', 'C=CCCCC'),
    ('2', 'CC(O)COC(C)COCc1ccc(cc1)C(OCC(C)OCC(C)O)CC'),
    ('3', 'FC(Cl)C(C)(C)OOC(=O)c1ccccc1'),
])
def test_gaff(name_smiles):
    name, smiles = name_smiles
    ff = ForceField.open(cwd + '/files/gaff.zfp')
    typer = ZftTyper(cwd + '/files/gaff.zft')
    mol = Molecule.from_smiles(f'{smiles} {name}')
    typer.type(mol)
    for atom in mol.atoms:
        print(atom)

    top = Topology([mol])
    system = System(top, ff, ignore_missing_improper=True)
    system.decompose_energy()
    system.minimize_energy()
