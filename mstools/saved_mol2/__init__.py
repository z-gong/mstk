import os

_PWD = os.path.dirname(os.path.abspath(__file__))

smiles_mol2_dict = {
    'F[P-](F)(F)(F)(F)F': os.path.join(_PWD, 'PF6-.mol2'),
    'F[P-](C(C(F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F': os.path.join(_PWD, 'PF3_C2F5_3-.mol2'),
    'CNC=O': os.path.join(_PWD, '03-C2H5NO.mol2'),
    'CNC(=O)C': os.path.join(_PWD, '04-C3H7NO.mol2'),
    'CCNC=O': os.path.join(_PWD, '05-C3H7NO.mol2'),
}
