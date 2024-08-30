import os
import tempfile
import filecmp
import pytest

from mstk.topology import Molecule, Topology
from mstk.forcefield import ForceField, Zfp, ZftTyper, GaffTyper
from mstk.simsys import System

cwd = os.path.dirname(os.path.abspath(__file__))

name_smiles_types = [
    ('a-01', 'CCCF', 'c3-c3-c3-f-hc-hc-hc-hc-hc-h1-h1'),
    ('a-02', 'ClC=CBr', 'cl-c2-c2-br-h1-h1'),
    ('a-03', 'COc1ccc(O)cc1', 'c3-os-ca-ca-ca-ca-oh-ca-ca-h1-h1-h1-ha-ha-ho-ha-ha'),
    ('a-04', 'ICC=O', 'i-c3-c-o-hc-hc-h1'),  # antechamber type the last hydrogen as ha
    ('a-05', 'CC(=O)CS', 'c3-c-o-c3-sh-hc-hc-hc-hc-hc-hs'),
    ('a-06', 'CC(=S)O', 'c3-c-s2-oh-hc-hc-hc-ho'),
    ('a-07', 'CSCC#C', 'c3-ss-c3-c1-c1-hc-hc-hc-hc-hc-hc'),
    ('a-08', 'CSSC', 'c3-ss-ss-c3-hc-hc-hc-hc-hc-hc'),
    ('a-09', 'CS(=O)C', 'c3-s4-o-c3-hc-hc-hc-hc-hc-hc'),
    ('a-10', 'CS(=O)(=O)NC', 'c3-s6-o-o-n-c3-hc-hc-hc-hn-h1-h1-h1'),
    ('a-11', 'CC(=O)N', 'c3-c-o-n-hc-hc-hc-hn-hn'),
    ('a-12', 'NCC=NO', 'n3-c3-c2-n2-oh-hn-hn-h1-h1-h1-ho'),
    ('a-13', 'C[NH3+]', 'c3-n4-h1-h1-h1-hn-hn-hn'),
    ('a-14', 'Cn1cccc1', 'c3-nb-ca-ca-ca-ca-h1-h1-h1-h4-ha-ha-h4'),  # antechamber does not consider pyrrole as aromatic
    ('a-15', 'c1ccccc1N', 'ca-ca-ca-ca-ca-ca-nh-ha-ha-ha-ha-ha-hn-hn'),
    ('a-16', 'CP=CCN(=O)=O', 'c3-p2-c2-c3-no-o-o-hc-hc-hc-hc-h1-h1'),
    ('a-17', 'CP', 'c3-p3-hc-hc-hc-hp-hp'),
    ('a-18', 'C[P](=O)C', 'c3-p4-o-c3-hc-hc-hc-hc-hc-hc'),
    ('a-19', 'CP(=O)(O)OC', 'c3-p5-o-oh-os-c3-hc-hc-hc-ho-h1-h1-h1'),
    ('b-01', 'C=CC=C', 'c2-ce-ce-c2-hc-hc-hc-hc-hc-hc'),
    ('b-02', 'C=CN=C', 'c2-ce-ne-c2-hc-hc-h1-h1-h1'),
    ('b-03', 'C=CC=CC=CC=C', 'c2-ce-ce-cf-cf-ce-ce-c2-hc-hc-hc-hc-hc-hc-hc-hc-hc-hc'),
    ('b-04', 'C#CC=C', 'c1-cg-ce-c2-hc-hc-hc-hc'),
    ('b-05', 'C=CP=C', 'c2-ce-pe-c2-hc-hc-hc-hc-hc'),
    ('b-06', 'o1cccc1', 'os-ca-ca-ca-ca-h4-ha-ha-h4'),  # antechamber does not consider furan as aromatic
    ('b-07', 'S1-C=C-N=C1', 'ss-ca-ca-nb-ca-ha-h4-h4'),  # antechamber does not consider 5-member ring as aromatic
    ('b-08', 'O1-C=P-C=C1', 'os-ca-pb-ca-ca-h4-ha-h4'),  # antechamber does not consider 5-member ring as aromatic
    ('b-09', 'C1=CNC=CC1=O', 'ca-ca-nb-ca-ca-ca-o-ha-h4-hn-h4-ha'),  # antechamber does not consider it as aromatic
    ('b-10', 'c1ccccc1c2ccccc2', 'ca-ca-ca-ca-ca-cp-cp-ca-ca-ca-ca-ca-ha-ha-ha-ha-ha-ha-ha-ha-ha-ha'),
    ('b-11', 'c1ccccc1[n+]2ccccc2', 'ca-ca-ca-ca-ca-cp-nb-ca-ca-ca-ca-ca-ha-ha-ha-ha-ha-h4-ha-ha-ha-h4'),
    ('b-12', 'c1ccccc1c2ccccc2c3ccccc3',
     'ca-ca-ca-ca-ca-cp-cp-ca-ca-ca-ca-cq-cq-ca-ca-ca-ca-ca-ha-ha-ha-ha-ha-ha-ha-ha-ha-ha-ha-ha-ha-ha'),
    ('b-13', 'c1ccccc1c2ccccc2c3cccc(c4ccco4)c3',
     'ca-ca-ca-ca-ca-cp-cp-ca-ca-ca-ca-cq-cq-ca-ca-ca-cp-cp-ca-ca-ca-os-ca-ha-ha-ha-ha-ha-ha-ha-ha-ha-ha-ha-ha-ha-ha-h4-ha'),
    # antechamer does not consider furan as aromatic
    ('b-14', 'c1cccc(cccc2)c12', 'ca-ca-ca-ca-ca-ca-ca-ca-ca-ca-ha-ha-ha-ha-ha-ha-ha-ha'),
    ('b-15', 'c1c(cccc3)c3cc(cccc2)c12', 'ca-ca-ca-ca-ca-ca-ca-ca-ca-ca-ca-ca-ca-ca-ha-ha-ha-ha-ha-ha-ha-ha-ha-ha'),
    ('b-16', 'CS(=C)C', 'c3-s4-c2-c3-hc-hc-hc-hc-hc-hc-hc-hc'),
    ('b-17', 'CS(=O)C=C', 'c3-sx-o-ce-c2-hc-hc-hc-hc-hc-hc'),
    ('b-18', 'C1CC1S(=O)(=O)C2CC=C2', 'cx-cx-cx-s6-o-o-cy-cy-cv-cv-hc-hc-hc-hc-hc-hc-hc-hc-hc-hc'),
    ('b-19', 'C1=CC1S(=O)(=O)C=C', 'cu-cu-cx-sy-o-o-ce-c2-hc-hc-hc-hc-hc-hc'),
    ('b-20', 'C[P](=S)O', 'c3-p4-s2-oh-hc-hc-hc-ho'),
    ('b-21', 'C[P](=O)C=C', 'c3-px-o-ce-c2-hc-hc-hc-hc-hc-hc'),
    ('b-22', 'CP(=N)(OC)OC', 'c3-p5-n2-os-c3-os-c3-hc-hc-hc-hn-h1-h1-h1-h1-h1-h1'),
    ('b-23', 'COP(=N)(OC)N=C', 'c3-os-py-n2-os-c3-ne-c2-h1-h1-h1-hn-h1-h1-h1-h1-h1'),
    ('b-24', 'COP(=O)(OC)P(=O)(OC)OC', 'c3-os-py-o-os-c3-py-o-os-c3-os-c3-h1-h1-h1-h1-h1-h1-h1-h1-h1-h1-h1-h1'),
    ('b-25', 'COP(O)(OC)=P(O)(OC)OC', 'c3-os-p5-oh-os-c3-p5-oh-os-c3-os-c3-h1-h1-h1-ho-h1-h1-h1-ho-h1-h1-h1-h1-h1-h1'),
]


@pytest.mark.parametrize('name_smiles_types', name_smiles_types)
def test_gaff_typer(name_smiles_types):
    name, smiles, str_types = name_smiles_types
    typer = GaffTyper()
    mol = Molecule.from_smiles(f'{smiles} {name}')
    typer.type(mol)

    assert '-'.join(atom.type for atom in mol.atoms) == str_types


@pytest.mark.skip(reason='many failures as expected, because of parameters missing in GAFF')
@pytest.mark.parametrize('name_smiles_types', name_smiles_types)
def test_gaff(name_smiles_types):
    name, smiles, str_types = name_smiles_types
    typer = GaffTyper()
    mol = Molecule.from_smiles(f'{smiles}')
    typer.type(mol)

    ff = ForceField.open('gaff.zff')
    top = Topology([mol])
    system = System(top, ff)
    system.decompose_energy()
    system.minimize_energy()


@pytest.mark.parametrize('smiles', [
    'CCC',
    'CCO',
    'c1ccccc1N',
    'c1ccccc1C(=O)N',
    'O=C=NCCCCCCN=C=O',
    'O=C=Nc1ccc(N=C=O)cc1',
])
def test_gaff_mod(smiles):
    typer = GaffTyper()
    mol = Molecule.from_smiles(f'{smiles}')
    typer.type(mol)

    ff = ForceField.open('gaff.zff')
    top = Topology([mol])
    system = System(top, ff)
    system.decompose_energy()
    system.minimize_energy()
