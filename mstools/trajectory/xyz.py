import numpy as np
from ..topology import Topology
from . import Trajectory, Frame


class Xyz(Trajectory):
    '''
    Since xyz format is not very useful, I only parse the first frame
    '''

    def __init__(self, file, mode='r'):
        super().__init__()
        if mode not in ('r', 'w', 'a'):
            raise Exception('Invalid mode')
        self._file = open(file, mode)
        if mode == 'r' or mode == 'a':
            self._get_info()

    def _get_info(self):
        '''
        Read the number of atoms and record the offset of lines and frames,
        so that we can read arbitrary frame later
        '''
        try:
            self.n_atom = int(self._file.readline())
        except:
            print('Invalid XYZ file')
            raise
        self._file.seek(0)

        self.n_frame = 1
        self._frame = Frame(self.n_atom)

    def read_frame(self, i_frame):
        if i_frame != 0:
            raise Exception('Can only read the first frame')

        return self._read_frame_from_string(self._file.read())

    def read_frames(self, i_frames: [int]):
        raise Exception('Can only read the first frame, use read_frame(0)')

    def _read_frame_from_string(self, string: str):
        frame = self._frame
        lines = string.splitlines()

        for i in range(self.n_atom):
            words = lines[i + 2].split()
            frame.positions[i] = np.array(list(map(float, words[1:4]))) / 10  # convert from A to nm

        return frame

    def write_frame(self, frame: Frame, topology: Topology, subset=None):
        if subset is None:
            subset = list(range(topology.n_atom))

        self._file.write('%i\n' % len(subset))
        self._file.write('Created by mstools: step= %i, t= %f ps\n' % (frame.step, frame.time))

        for ii, id in enumerate(subset):
            atom = topology.atoms[id]
            pos = frame.positions[id] * 10  # convert from nm to A
            # atom.symbol is more friendly than atom.type for visualizing trajectory
            line = '%-8s %10.5f %10.5f %10.5f\n' % (atom.symbol, pos[0], pos[1], pos[2])
            self._file.write(line)
