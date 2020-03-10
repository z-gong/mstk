import numpy as np
from ..topology import Topology
from . import Trajectory, Frame


class Dcd(Trajectory):
    '''
    Currently mstools use mdtraj to parse DCD format
    '''

    def __init__(self, file, mode='r'):
        super().__init__()

        try:
            from mdtraj.formats import DCDTrajectoryFile
        except:
            raise ImportError('Currently mstools use mdtraj to parse DCD format. Cannot import mdtraj')

        self._dcd = DCDTrajectoryFile(file)

        if mode == 'r':
            self._get_info()
        elif mode == 'w':
            pass

    def close(self):
        self._dcd.close()

    def _get_info(self):
        '''
        Read the number of atoms and record the offset of lines and frames,
        so that we can read arbitrary frame later
        '''
        self._mdtraj_positions, self._mdtraj_cell_lengths, self._mdtraj_cell_angles = self._dcd.read()
        self.n_frame, self.n_atom, _ = self._mdtraj_positions.shape
        if np.any(self._mdtraj_cell_angles != 90):
            raise Exception('Currently non-rectangular box is not supported')

    def read_frame(self, i_frame):
        return self.read_frames([i_frame])[0]

    def read_frames(self, i_frames: [int]) -> [Frame]:
        frames = []
        for i in i_frames:
            frame = Frame(self.n_atom)
            frame.positions = self._mdtraj_positions[i] / 10  # convert A to nm
            frame.box = self._mdtraj_cell_lengths[i] / 10  # convert A to nm
            frames.append(frame)
        return frames

    def write_frame(self, topology, frame, subset=None):
        if subset is None:
            subset = list(range(topology.n_atom))
        self._dcd.write(frame.positions[subset] * 10, frame.box * 10, np.array([90, 90, 90]))  # convert nm to A
