import numpy as np
from . import Trajectory, Frame


class Xtc(Trajectory):
    '''
    Currently mstools use mdtraj to parse XTC format
    '''

    def __init__(self, file, mode='r'):
        super().__init__()

        try:
            from mdtraj.formats import XTCTrajectoryFile
        except:
            raise ImportError(
                'Currently mstools use mdtraj to parse XTC format. Cannot import mdtraj')

        if mode == 'r':
            self._xtc = XTCTrajectoryFile(file, mode='r')
            self._get_info()
        elif mode == 'w':
            self._xtc = XTCTrajectoryFile(file, mode='w')
        else:
            raise Exception('Appending not supported for XTC')

    def close(self):
        try:
            self._xtc.close()
        except:
            pass

    def _get_info(self):
        '''
        Read the number of atoms and record the offset of lines and frames,
        so that we can read arbitrary frame later
        '''
        self.n_frame = len(self._xtc)
        if self.n_frame == 0:
            raise Exception('Empty XTC file')
        positions, times, steps, box_vectors = self._xtc.read(1)
        _, self.n_atom, _ = positions.shape
        self._frame = Frame(self.n_atom)

    def read_frame(self, i_frame):
        if i_frame >= self.n_frame:
            raise Exception('i_frame should be smaller than %i' % self.n_frame)
        frame = self._frame
        self._xtc.seek(i_frame)
        positions, times, steps, box_vectors = self._xtc.read(1)
        vector = box_vectors[0]
        vector[np.abs(vector) < 1E-4] = 0  # in case precision issue
        frame.positions = positions[0]
        frame.cell.set_box(vector)
        frame.time = times[0]
        frame.step = steps[0]

        return frame

    def read_frames(self, i_frames: [int]) -> [Frame]:
        if any(i >= self.n_frame for i in i_frames):
            raise Exception('i_frame should be smaller than %i' % self.n_frame)
        frames = []
        for i in i_frames:
            frame = Frame(self.n_atom)
            self._xtc.seek(i)
            positions, times, steps, box_vectors = self._xtc.read(1)
            box_vector = box_vectors[0]
            box_vector[np.abs(box_vector) < 1E-4] = 0  # in case precision issue
            frame.positions = positions[0]
            frame.cell.set_box(box_vector)
            frames.append(frame)

        return frames

    def write_frame(self, frame, topology=None, subset=None):
        if subset is None:
            positions = frame.positions
        else:
            positions = frame.positions[subset]
        self._xtc.write(positions, frame.time, frame.step, frame.cell.vectors)
