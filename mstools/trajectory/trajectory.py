import numpy as np


class Frame():
    def __init__(self, n_atom):
        self.step = 0
        self._box = np.array([0., 0., 0.])
        self.positions = np.array([[0., 0., 0.]] * n_atom)
        self.has_velocity = False
        self.velocities = np.array([[0., 0., 0.]] * n_atom)
        self.has_charge = False
        self.charges = np.array([0.] * n_atom)  # for fluctuating charge simulations

    @property
    def box(self):
        return self._box

    @box.setter
    def box(self, value):
        if not isinstance(value, (list, tuple, np.ndarray)) or len(value) != 3:
            raise ValueError('box should has three elements')
        self._box = np.array(value)


class Trajectory():
    def __init__(self):
        self.n_atom = 0
        self.n_frame = 0

    @staticmethod
    def open(file, mode='r'):
        from .gro import Gro
        from .pdb import Pdb
        from .xyz import Xyz
        from .lammps import LammpsTrj

        if file.endswith('.lammpstrj') or file.endswith('.ltrj'):
            return LammpsTrj(file, mode)
        elif file.endswith('.gro'):
            return Gro(file, mode)
        elif file.endswith('.pdb'):
            return Pdb(file, mode)
        elif file.endswith('.xyz'):
            return Xyz(file, mode)
        else:
            raise Exception('filename for trajectory not understand')
