import numpy as np

from .topology import Atom, Molecule


class Frame():
    def __init__(self, n_atom):
        self.step = 0
        self.xlo = self.xhi = self.ylo = self.yhi = self.zlo = self.zhi = 0.0
        self.box = np.array([0., 0., 0.])
        self.positions = np.array([[0., 0., 0.]] * n_atom)
        self.charges = np.array([0.] * n_atom)  # for fluctuating charge simulations


class Trajectory():
    def __init__(self):
        self.n_atom = 0
        self.n_frame = 0
