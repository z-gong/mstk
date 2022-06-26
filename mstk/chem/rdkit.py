from mstk.errors import RDKitError

try:
    from rdkit.Chem import AllChem as Chem
except ImportError:
    RDKIT_FOUND = False
else:
    RDKIT_FOUND = True


def create_mol_from_smiles(smiles):
    '''
    Create a RDKit molecule object from SMILES string.

    Parameters
    ----------
    smiles : str

    Returns
    -------
    rdmol : rdkit.rdchem.Mol
    '''
    if not RDKIT_FOUND:
        raise ImportError('RDKit is required for parsing SMILES')

    try:
        rdmol = Chem.MolFromSmiles(smiles)
    except:
        raise RDKitError('Invalid SMILES')

    rdmol = Chem.AddHs(rdmol)
    if Chem.EmbedMolecule(rdmol, useRandomCoords=True) == -1:
        if Chem.EmbedMolecule(rdmol) == -1:
            raise RDKitError('Failed generating 3D coordinates')

    return rdmol
