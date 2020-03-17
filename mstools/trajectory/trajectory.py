from io import IOBase
import numpy as np
from ..topology import Topology, UnitCell


class Frame():
    def __init__(self, n_atom):
        self.step = 0
        self.unitcell: UnitCell = None
        self.positions = np.zeros((n_atom, 3), dtype=np.float32)
        self.has_velocity = False
        self.velocities = np.zeros((n_atom, 3), dtype=np.float32)
        self.has_charge = False
        self.charges = np.zeros(n_atom, dtype=np.float32)  # for fluctuating charge simulations

    @property
    def box(self):
        if self.unitcell is None:
            return None
        return self.unitcell.box

    @box.setter
    def box(self, value):
        if self.unitcell is not None and not self.unitcell.is_rectangular:
            raise Exception('unitcell is not rectangular, set unitcell vectors instead of box')

        if not isinstance(value, (list, tuple, np.ndarray)) or len(value) != 3:
            raise ValueError('box should has three elements')

        if self.unitcell is None:
            self.unitcell = UnitCell(value)
        else:
            self.unitcell.vectors = value


class Trajectory():
    def __init__(self):
        self.n_atom = 0
        self.n_frame = 0
        self._file = IOBase()
        self._frame: Frame = None

    def __del__(self):
        self.close()

    def close(self):
        self._file.close()

    def read_frame(self, i_frame: int):
        '''
        Read a single frame
        @param i_frame:
        @return:
        '''
        pass

    def read_frames(self, i_frames: [int]):
        '''
        Read a bunch of frames for multiprocessing
        @param i_frames:
        @return:
        '''
        pass

    def write_frame(self, topology: Topology, frame: Frame, subset: [int]):
        '''
        Write one frame to trajectory file
        # TODO fix the performance of subset
        '''
        pass

    @staticmethod
    def open(file, mode='r'):
        from .gro import Gro
        from .pdb import Pdb
        from .xyz import Xyz
        from .lammps import LammpsTrj
        from .dcd import Dcd
        from .combined_trajectory import CombinedTrajectory

        if isinstance(file, list):
            return CombinedTrajectory(file, mode)
        elif file.endswith('.lammpstrj') or file.endswith('.ltrj'):
            return LammpsTrj(file, mode)
        elif file.endswith('.gro'):
            return Gro(file, mode)
        elif file.endswith('.pdb'):
            return Pdb(file, mode)
        elif file.endswith('.xyz'):
            return Xyz(file, mode)
        elif file.endswith('.dcd'):
            return Dcd(file, mode)
        else:
            raise Exception('filename for trajectory not understand')
