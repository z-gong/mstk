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
        # assume all frames have the same number of atoms, but n_atom can only be determined later
        self.n_atom = -1  # -1 means unknown yet
        self.n_frame = -1  # -1 means unknown yet
        self._file = IOBase()
        self._mode = 'r'
        self._opened = False
        # we don't have to create a new Frame object every time we call read_frame()
        self._frame: Frame = None

    def __del__(self):
        self.close()

    def close(self):
        self._file.close()
        self._opened = False

    def read_frame(self, i_frame: int):
        '''
        Read a single frame
        @param i_frame:
        @return:
        '''
        if self._mode != 'r' or not self._opened:
            raise Exception('mode != "r" or closed trajectory')

        if i_frame >= self.n_frame:
            raise Exception('i_frame should be smaller than %i' % self.n_frame)

        if self._frame is None:
            if self.n_atom == -1:
                raise Exception('Invalid number of atoms')
            self._frame = Frame(self.n_atom)

        self._read_frame(i_frame, self._frame)
        return self._frame

    def read_frames(self, i_frames: [int]):
        '''
        Read a bunch of frames for multiprocessing
        @param i_frames:
        @return:
        '''
        if self._mode != 'r' or not self._opened:
            raise Exception('mode != "r" or closed trajectory')

        if any(i >= self.n_frame for i in i_frames):
            raise Exception('i_frame should be smaller than %i' % self.n_frame)

        frames = [Frame(self.n_atom) for _ in i_frames]
        for ii, i_frame in enumerate(i_frames):
            self._read_frame(i_frame, frames[ii])
        return frames

    def write_frame(self, frame: Frame, topology: Topology = None, subset: [int] = None, **kwargs):
        '''
        Write one frame to trajectory file
        '''
        if self._mode not in ('w', 'a') or not self._opened:
            raise Exception('mode not in ("w", "a") or closed trajectory')

        self._write_frame(frame, topology, subset, **kwargs)

    def _read_frame(self, i_frame, frame):
        raise NotImplementedError('Method not implemented')

    def _write_frame(self, frame, topology, subset, **kwargs):
        raise NotImplementedError('Method not implemented')

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
