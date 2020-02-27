import numpy as np


class Atom():
    def __init__(self, name='UNK'):
        self.id = -1
        self.molecule: Molecule = None
        self.name = name
        self.type = ''
        self.symbol = ''
        self.charge = 0.
        self.mass = 0.
        self.position = np.array([0, 0, 0], dtype=float)

    def __repr__(self):
        return f'<Atom: {self.name} {self.id} {self.type}>'


class Molecule():
    def __init__(self, name='UNK'):
        self.id = -1
        self.name = name
        self.atoms: [Atom] = []

    def add_atom(self, atom: Atom):
        self.atoms.append(atom)
        atom.molecule = self

    def remove_atom(self, atom: Atom):
        self.atoms.remove(atom)
        atom.molecule = None

    def __repr__(self):
        return f'<Molecule: {self.name} {self.id}>'


class Topology():
    def __init__(self):
        self.n_atom = 0
        self.atoms: [Atom] = []
        self.n_molecule = 0
        self.molecules: [Molecule] = []
        self.remark = ''
        self.is_drude = False

    def init_from_molecules(self, molecules: [Molecule]):
        self.n_molecule = len(molecules)
        self.molecules: [Molecule] = molecules[:]
        self.assign_id()

    def assign_id(self):
        idx_atom = 0
        for i, mol in enumerate(self.molecules):
            mol.id = i
            for j, atom in mol.atoms:
                atom.id = idx_atom
                idx_atom += 1
