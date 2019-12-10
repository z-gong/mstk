import numpy as np


class Atom():
    def __init__(self, id=0, name='UNK'):
        self.id = id
        self.molecule: Molecule = None
        self.name = name
        self.type = ''
        self.element = ''
        self.charge = 0.
        self.position = np.array([0, 0, 0], dtype=float)


class Molecule():
    def __init__(self, id=0, name='UNK'):
        self.id = id
        self.name = name
        self.atoms: [Atom] = []

    def add_atom(self, atom: Atom):
        self.atoms.append(atom)
        atom.molecule = self

    def remove_atom(self, atom: Atom):
        self.atoms.remove(atom)
        atom.molecule = None
