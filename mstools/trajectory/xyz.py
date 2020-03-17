import numpy as np
from ..topology import Topology
from . import Trajectory, Frame


class Xyz(Trajectory):
    '''
    Since xyz format is not very useful, I only parse the first frame
    '''

    def __init__(self, file, mode='r'):
        super().__init__()
        self._file = open(file, mode)
        if mode == 'r':
            self._get_info()
        elif mode == 'w':
            pass

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

    def write_frame(self, topology: Topology, frame: Frame, subset=None):
        self._file.write('%i\n' % topology.n_atom)
        self._file.write('\n')

        if subset is None:
            subset = list(range(topology.n_atom))
        for ii, id in enumerate(subset):
            atom = topology.atoms[id]
            pos = frame.positions[id] * 10  # convert from nm to A
            line = '%-8s %10.5f %10.5f %10.5f\n' % (atom.type, pos[0], pos[1], pos[2])
            self._file.write(line)
