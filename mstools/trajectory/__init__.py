import numpy as np

from ..molio import Atom


class Frame():
    def __init__(self, n_atom):
        self.step = 0
        self.atoms = np.array([Atom() for i in range(n_atom)], dtype=object)

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
        self.current_frame: Frame = None

    def read_next_frame(self):
        pass
