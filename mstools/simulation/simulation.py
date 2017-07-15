import math

from ..utils import create_mol_from_smiles, estimate_density_from_formula
from ..wrapper import Packmol, DFF


class Simulation():
    def __init__(self, packmol_bin=None, dff_root=None, jobmanager=None):
        self.packmol = Packmol(packmol_bin=packmol_bin)
        self.dff = DFF(dff_root=dff_root)
        self.jobmanager = jobmanager
        self.procedure = None

        self.n_mol_list: [int]
        self.msd: str = 'init.msd'

    def set_procedure(self, procedure):
        self.procedure = procedure

    def build(self):
        pass

    def prepare(self):
        pass

    def run(self):
        self.jobmanager.submit()

    def check_finished(self):
        pass

    def analyze(self) -> {str: float}:
        pass

    def set_system(self, smiles_list: [str], n_atoms: int, density: float = None):
        self.pdb_list = []
        self.mol2_list = []
        n_components = len(smiles_list)
        n_atom_list = []  # number of atoms of each molecule
        molwt_list = []  # molecule weight of each molecule
        density_list = []  # estimated density of each molecule
        for i, smiles in enumerate(smiles_list):
            pdb = 'mol-%i.pdb' % (i + 1)
            mol2 = 'mol-%i.mol2' % (i + 1)
            py_mol = create_mol_from_smiles(smiles, pdb_out=pdb, mol2_out=mol2)
            self.pdb_list.append(pdb)
            self.mol2_list.append(mol2)
            n_atom_list.append(len(py_mol.atoms))
            molwt_list.append(py_mol.molwt)
            density_list.append(estimate_density_from_formula(py_mol.formula))

        self.n_mol_list = [math.ceil(n_atoms / n_components / n_atom) for n_atom in n_atom_list]
        mass = sum([molwt_list[i] * self.n_mol_list[i] for i in range(n_components)])
        if density is None:
            density = sum([density_list[i] * self.n_mol_list[i] for i in range(n_components)]) / sum(self.n_mol_list)
        self.length = (10 / 6.022 * mass / density) ** (1 / 3)  # assume cubic box
