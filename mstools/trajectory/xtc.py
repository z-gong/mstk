import numpy as np
from . import Trajectory, Frame

try:
    import chemfiles
except ImportError:
    _CHEMFILES_FOUND = False
else:
    _CHEMFILES_FOUND = True


class Xtc(Trajectory):
    '''
    Currently mstools use mdtraj to parse XTC format
    '''

    def __init__(self, file, mode='r'):
        super().__init__()

        if not _CHEMFILES_FOUND:
            raise ImportError(
                'Currently mstools use chemfiles to parse XTC format. Cannot import chemfiles')

        if mode not in ('r', 'w', 'a'):
            raise Exception('Invalid mode')

        self._xtc = chemfiles.Trajectory(file, mode)
        if mode == 'r':
            self._get_info()

        self._mode = mode
        self._opened = True

    def close(self):
        try:
            self._xtc.close()
        except:
            pass
        self._opened = False

    def _get_info(self):
        '''
        Read the number of atoms and record the offset of lines and frames,
        so that we can read arbitrary frame later
        '''
        self.n_frame = self._xtc.nsteps
        if self.n_frame == 0:
            raise Exception('Empty XTC file')
        cf_frame = self._xtc.read()
        self.n_atom = len(cf_frame.atoms)
        self._frame = Frame(self.n_atom)

    def _read_frame(self, i_frame, frame):
        cf_frame = self._xtc.read_step(i_frame)
        cf_cell = cf_frame.cell
        lengths = [i / 10 for i in cf_cell.lengths]
        angles = cf_cell.angles
        frame.positions = cf_frame.positions.astype(np.float32) / 10
        frame.cell.set_box([lengths, angles])
        frame.step = frame.step

    def _write_frame(self, frame, topology=None, subset=None, **kwargs):
        if subset is None:
            positions = frame.positions
        else:
            positions = frame.positions[subset]

        cf_frame = chemfiles.Frame()
        cf_frame.resize(len(positions))
        cf_frame.positions[:] = positions * 10
        cf_frame.cell = chemfiles.UnitCell(*(frame.cell.lengths * 10), *frame.cell.angles)
        cf_frame.step = frame.step

        self._xtc.write(cf_frame)
