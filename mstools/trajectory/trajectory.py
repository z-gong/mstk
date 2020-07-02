import numpy as np
from io import IOBase
from ..topology import Topology, UnitCell
from .. import logger


class Frame():
    '''
    A frame records the positions and unit cell (and optionally velocities, charges, current step and current time)
    of the system at one specific time step.

    Velocity is stored for analyzing dynamic properties with auto correlation function.
    Charge is stored for fluctuating charge method.

    Parameters
    ----------
    n_atom : int
        The number of atoms in the frame.

    Attributes
    ----------
    step : int
        The current step. -1 means unknown.
    time : float
        The current simulation time. -1 means unknown.
    cell : UnitCell
        The unit cell.
    positions : array_like
        The positions of all atoms in the frame.
        It's a np.ndarry of shape (n_atom, 3)
    has_velocity : bool
        Whether or not there's velocity information for all atoms
    velocities : array_like
        The velocities of all atoms in the frame.
        It's a np.ndarry of shape (n_atom, 3)
    has_charge : bool
        Whether or not there's charge information for all atoms
    charges : array_like
        The charge of all atoms in the frame.
        It's a np.ndarry of shape (n_atom,)
    '''
    def __init__(self, n_atom):
        self.step = -1  # -1 means step information not found
        self.time = -1  # in ps. -1 means time information not found
        self.cell = UnitCell([0., 0., 0.])
        self.positions = np.zeros((n_atom, 3), dtype=np.float32)
        self.has_velocity = False
        self.velocities = np.zeros((n_atom, 3), dtype=np.float32)
        self.has_charge = False
        self.charges = np.zeros(n_atom, dtype=np.float32)  # for fluctuating charge simulations

    def resize(self, n_atom):
        '''
        Resize the number of atoms recorded in the frame.

        Usually it is not necessary to call this method explicitly.
        If n_atom is smaller than previous value, the information of the end atoms will be trimmed.
        If n_atom is larger than previous value, the information for the new atoms will be random .

        Parameters
        ----------
        n_atom : int
            The new number of atoms
        '''
        if n_atom < len(self.positions):
            logger.warning('n_atom is smaller than original. '
                           'The positions, velocities and charges at the end will be lost')
        self.positions.resize((n_atom, 3), refcheck=False)
        self.velocities.resize((n_atom, 3), refcheck=False)
        self.charges.resize((n_atom, 3), refcheck=False)


class Trajectory():
    '''
    Base class of Trajectory reader/writers.

    A Trajectory is made of a series of Frames.
    It is assumed that all frames have the same number of atoms.
    If it's not the case (e.g. the LammpsDump can have variable number of atoms in each frame),
    then this class may not be able to work correctly.

    The frames are loaded on request for performance issue,
    therefore a file handler is kept open until :func:`close` is called.

    The Trajectory should not be constructed with its :func:`Trajectory.__init__` function.
    Instead, use the :func:`open` static function to open a trajectory.
    A suitable parser will be determined automatically from the extension of the file.
    Or, explicitly call a subclass of Trajectory to open a file.

    Several files can be opened at the same time (even for files of different format),
    which is useful for processing truncated trajectories.
    In this case, the trajectory is read-only.

    Attributes
    ----------
    n_atom : int
        Number of atoms in the trajectory. -1 means unknown yet
    n_frame : int
        Number of frames in the trajectory. -1 means unknown yet

    Examples
    --------
    >>> trj = Trajectory.open('input.dcd')
    >>> frame = trj.read_frame(0)
    >>> trj.close()

    >>> trj = Trajectory.open('output.dcd', mode='w')
    >>> trj.write(frame)
    >>> trj.close()

    >>> trj = Trajectory.open('output.dcd', mode='a')
    >>> trj.write(frame)
    >>> trj.close()

    >>> trj = Trajectory.open(['input1.dcd', 'input2.xtc'])
    '''

    def __init__(self):
        self.n_atom = -1
        self.n_frame = -1
        self._file = IOBase()
        self._mode = 'r'
        self._opened = False
        # we don't have to create a new Frame object every time we call read_frame()
        self._frame: Frame = None

    def __del__(self):
        self.close()

    def close(self):
        '''
        Close the opened trajectory file(s).
        '''
        self._file.close()
        self._opened = False

    def read_frame(self, i_frame: int):
        '''
        Read a single frame from the trajectory.

        For the best of performance, a cached frame is used.
        When this method is called, the position and cell etc in the cached frame are updated,
        and the reference of this frame will be returned.
        Therefore, if you want to handle several different frames at the same time,
        this method should not be used because all the frames actually point to the same frame.
        In this case, use :func:`read_frames` instead.

        i_frame should be in the range of (0, n_frame), otherwise and Exception will be raised.

        Parameters
        ----------
        i_frame : int
            The index of the frame in trajectory.

        Returns
        -------
        frame : Frame

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
        Read a bunch of frames from the trajectory.

        Unlike :func:`read_frame`, the cached frame is not used.
        Instead, a new Frame object is constructed for each frame.
        This method should be called when you want to multiprocess several frames in parallel.

        All the items in i_frames should be in the range of (0, n_frame), otherwise and Exception will be raised.

        Parameters
        ----------
        i_frames : list of int

        Returns
        -------
        frames : list of Frame

        '''
        if self._mode != 'r' or not self._opened:
            raise Exception('mode != "r" or closed trajectory')

        if any(i >= self.n_frame for i in i_frames):
            raise Exception('i_frame should be smaller than %i' % self.n_frame)

        frames = [Frame(self.n_atom) for _ in i_frames]
        for ii, i_frame in enumerate(i_frames):
            self._read_frame(i_frame, frames[ii])
        return frames

    def write_frame(self, frame, topology=None, subset=None, **kwargs):
        '''
        Write one frame to the opened trajectory file.

        It requires the trajectory opened in 'w' or 'a' mode.
        Some trajectory format (e.g. GRO, XYZ) record part of the topology information,
        therefore, a topology should be provided for these formats.
        While it's not required for more space efficient formats like DCD, XTC.

        It is possible to write only a subset of atoms in the frame to the trajectory file.
        In this case, make sure the number of atoms in the subset is the same for every frame been written.

        Parameters
        ----------
        frame : Frame
        topology : Topology, optional
        subset : list of int, optional
            The index of atoms in the frame to be written
        kwargs : dict
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
        '''
        Open (a) trajectory file(s)

        The file(s) can be opened in 'r', 'w' or 'a' modes, which mean read, write and append.
        The class used for parsing the trajectory is determined by the extension of the file.

        Parameters
        ----------
        file : str or list of str
            If it's a str, then open the file with corresponding parser.
            if it's a list of str, then call :class:`CombinedTrajectory` to open all of the files.
        mode : ['r', 'w', 'a']

        Returns
        -------
        trajectory : subclass of Trajectory

        '''
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

    @staticmethod
    def read_frame_from_file(file, i_frame=0):
        '''
        Read one single frame from a trajectory file.

        This is simply a shortcut for

        >>> trj = Trajectory.open('input.dcd')
        >>> frame = trj.read_frame(i_frame)
        >>> trj.close()

        By default, the first frame is read.
        Otherwise, argument `i_frame` should be set.
        If i_frame set to -1, then the last frame is read.

        Parameters
        ----------
        file : str
        i_frame : int

        Returns
        -------
        frame : Frame

        '''
        trj = Trajectory.open(file, 'r')
        if i_frame == -1:
            i_frame = trj.n_frame - 1
        frame = trj.read_frame(i_frame)
        trj.close()
        return frame
