import numpy as np


class Frame():
    def __init__(self, n_atom):
        self.step = 0
        self.xlo = self.xhi = self.ylo = self.yhi = self.zlo = self.zhi = 0.0
        self.box = np.array([0., 0., 0.])
        self.positions = np.array([[0., 0., 0.]] * n_atom)
        self.has_charge = False
        self.charges = np.array([0.] * n_atom)  # for fluctuating charge simulations
        self.has_velocity = False
        self.velocities = np.array([[0., 0., 0.]] * n_atom)


class Trajectory():
    def __init__(self):
        self.n_atom = 0
        self.n_frame = 0

    @staticmethod
    def open(file, mode='r'):
        from .gro import GRO
        from .pdb import PDB
        from .xyz import XYZ
        from .lammps import LammpsTrj

        if file.endswith('.lammpstrj') or file.endswith('.ltrj'):
            return LammpsTrj(file, mode)
        elif file.endswith('.gro'):
            return GRO(file, mode)
        elif file.endswith('.pdb'):
            return PDB(file, mode)
        elif file.endswith('.xyz'):
            return XYZ(file, mode)
        else:
            raise Exception('filename for trajectory not understand')
