import os

_PWD = os.path.dirname(os.path.abspath(__file__))

smiles_mol2_dict = {
    'F[P-](F)(F)(F)(F)F': os.path.join(_PWD, 'PF6-.mol2'),
    'CNC=O': os.path.join(_PWD, 'methylformamide.mol2'),
    'CNC(=O)C': os.path.join(_PWD, 'methylacetamide.mol2'),
    'CCNC=O': os.path.join(_PWD, 'ethylformamide.mol2'),
    'CC1=CCCC(=CCC1)C': os.path.join(_PWD, 'C10H16.mol2'),
}
