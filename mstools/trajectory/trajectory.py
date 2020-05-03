from io import IOBase
import numpy as np
from ..topology import Topology, UnitCell


class Frame():
    def __init__(self, n_atom):
        self.step = -1  # -1 means step information not found
        self.time = -1  # in ps. -1 means time information not found
        self.cell = UnitCell([0., 0., 0.])
        self.positions = np.zeros((n_atom, 3), dtype=np.float32)
        self.has_velocity = False
        self.velocities = np.zeros((n_atom, 3), dtype=np.float32)
        self.has_charge = False
        self.charges = np.zeros(n_atom, dtype=np.float32)  # for fluctuating charge simulations


class Trajectory():
    def __init__(self):
        # TODO what if the number of atoms is different in each step
        self.n_atom = 0
        self.n_frame = 0
        self._file = IOBase()
        # so that we don't need to create a new Frame object every time we call read_frame()
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

    def write_frame(self, frame: Frame, topology: Topology, subset: [int]):
        '''
        Write one frame to trajectory file
        '''
        pass

    @staticmethod
    def open(file, mode='r'):
        modes_allowed = ('r', 'w', 'a')
        if mode not in modes_allowed:
            raise Exception('mode should be one of %s' % str(modes_allowed))

        from . import Gro, Xyz, LammpsTrj, Dcd, Xtc, CombinedTrajectory

        if isinstance(file, list):
            return CombinedTrajectory(file, mode)
        elif file.endswith('.lammpstrj') or file.endswith('.ltrj'):
            return LammpsTrj(file, mode)
        elif file.endswith('.gro'):
            return Gro(file, mode)
        elif file.endswith('.xyz'):
            return Xyz(file, mode)
        elif file.endswith('.dcd'):
            return Dcd(file, mode)
        elif file.endswith('.xtc'):
            return Xtc(file, mode)
        else:
            raise Exception('filename for trajectory not understand')
