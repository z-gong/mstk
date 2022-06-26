import numpy as np
from mstk.trajectory import Frame
from mstk.trajectory.handler import TrjHandler


class Dcd(TrjHandler):
    '''
    Read and write cell and positions from/to DCD file.

    Velocities, step and time are ignored.

    Currently mstk use mdtraj to support DCD file.
    Therefore, appending is not supported.
    '''

    def __init__(self, file, mode='r'):
        super().__init__()

        try:
            from mdtraj.formats import DCDTrajectoryFile
        except:
            raise ImportError('Currently mstk use mdtraj to parse DCD format. Cannot import mdtraj')

        if mode == 'r':
            self._dcd = DCDTrajectoryFile(file, mode='r')
        elif mode == 'w':
            self._dcd = DCDTrajectoryFile(file, mode='w')
        else:
            raise Exception('Appending not supported for DCD')

    def close(self):
        try:
            self._dcd.close()
        except:
            pass

    def get_info(self):
        self.n_frame = len(self._dcd)
        if self.n_frame == 0:
            raise Exception('Empty DCD file')
        positions, lengths, angles = self._dcd.read(1)
        _, self.n_atom, _ = positions.shape

        return self.n_atom, self.n_frame

    def read_frame(self, i_frame, frame):
        self._dcd.seek(i_frame)
        positions, box_lengths, box_angles = self._dcd.read(1)
        angle = box_angles[0]
        angle[np.abs(angle - 90) < 1E-4] = 90  # in case precision issue
        frame.positions = positions[0] / 10  # convert A to nm
        frame.cell.set_box([box_lengths[0] / 10, angle])  # convert A to nm

    def write_frame(self, frame, subset=None, **kwargs):
        '''
        Write a frame into the opened DCD file

        Parameters
        ----------
        frame : Frame
        subset : list of int, optional
        kwargs : dict
            Ignored
        '''
        if subset is None:
            positions = frame.positions
        else:
            positions = frame.positions[subset]
        self._dcd.write(positions * 10, frame.cell.lengths * 10, frame.cell.angles)  # convert nm to A

TrjHandler.register_format('.dcd', Dcd)
