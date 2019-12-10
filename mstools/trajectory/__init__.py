import numpy as np

from ..molio import Atom, Molecule


class Frame():
    def __init__(self, n_atom, n_molecule=None):
        self.step = 0
        self.xlo = self.xhi = self.ylo = self.yhi = self.zlo = self.zhi = 0.0
        self.box = np.array([0, 0, 0], dtype=float)
        self.atoms = np.array([Atom(id=i) for i in range(n_atom)], dtype=object)
        if n_molecule is None:
            n_molecule = n_atom
        self.molecules = np.array([Molecule(id=i) for i in range(n_molecule)], dtype=object)

    @property
    def positions(self):
        return np.array([atom.position for atom in self.atoms])

    @property
    def charges(self):
        return np.array([atom.charge for atom in self.atoms])


class Trajectory():
    def __init__(self):
        self.n_atom = 0
        self.n_frame = 0
