import numpy as np
from mstk.chem.constant import *
from mstk.trajectory import Frame
from mstk.trajectory.handler import TrjHandler


class Dcd(TrjHandler):
    '''
    Read and write step, cell and positions from/to DCD file.

    Currently mstk use chemfiles to support DCD file.
    '''

    def __init__(self, file, mode='r'):
        super().__init__()

        try:
            import chemfiles
        except:
            raise ImportError('Currently mstk use chemfiles to parse DCD format. Cannot import chemfiles')

        if mode not in ('r', 'w', 'a'):
            raise Exception('Invalid mode')

        self._dcd = chemfiles.Trajectory(file, mode, format='DCD')

    def close(self):
        try:
            self._dcd.close()
        except:
            pass

    def get_info(self):
        self.n_frame = self._dcd.nsteps
        if self.n_frame == 0:
            raise Exception('Empty DCD file')
        cf_frame = self._dcd.read()
        self.n_atom = len(cf_frame.atoms)

        return self.n_atom, self.n_frame

    def read_frame(self, i_frame, frame):
        cf_frame = self._dcd.read_step(i_frame)
        cf_cell = cf_frame.cell
        lengths = [i / 10 for i in cf_cell.lengths]
        angles = [i * DEG2RAD for i in cf_cell.angles]
        frame.positions = cf_frame.positions.astype(float) / 10
        frame.cell.set_box([lengths, angles])
        frame.step = frame.step

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
        import chemfiles

        if subset is None:
            positions = frame.positions
        else:
            positions = frame.positions[subset]

        cf_frame = chemfiles.Frame()
        cf_frame.resize(len(positions))
        cf_frame.positions[:] = positions * 10
        cf_frame.cell = chemfiles.UnitCell(frame.cell.lengths * 10, frame.cell.angles * RAD2DEG)
        cf_frame.step = frame.step

        self._dcd.write(cf_frame)


TrjHandler.register_format('.dcd', Dcd)
