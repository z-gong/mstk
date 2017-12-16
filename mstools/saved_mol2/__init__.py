import os

_PWD = os.path.dirname(os.path.abspath(__file__))

smiles_mol2_dict = {
    'F[P-](F)(F)(F)(F)F': os.path.join(_PWD, 'PF6-.mol2'),
    'F[P-](C(C(F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F': os.path.join(_PWD, 'PF3_C2F5_3-.mol2')
}
