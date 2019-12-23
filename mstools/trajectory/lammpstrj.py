import traceback
import numpy as np
from io import StringIO
import pandas as pd
from . import Trajectory, Frame
from .topology import Atom, Molecule


class LammpsTrj(Trajectory):
    def __init__(self, trj_file, mode='r'):
        super().__init__()
        self._file = open(trj_file, mode)
        self._get_info()

    def _get_info(self):
        '''
        Read the number of atoms and record the offset of lines and frames,
        so that we can read arbitrary frame later
        '''
        try:
            self._file.readline()
            self._file.readline()
            self._file.readline()
            self.n_atom = int(self._file.readline())
        except:
            print('Invalid lammpstrj file')
            raise
        self._file.seek(0)

        # read in the file once and build a list of line offsets
        self._line_offset = []
        offset = 0
        for line in self._file:
            self._line_offset.append(offset)
            offset += len(line)
        # the last element is the length of whole file
        self._line_offset.append(offset)
        self._file.seek(0)

        # build a list of frame offsets
        self.n_frame = len(self._line_offset) // (9 + self.n_atom)
        self._frame_offset = []
        for i in range(self.n_frame + 1):
            line_start = (9 + self.n_atom) * i
            self._frame_offset.append(self._line_offset[line_start])

    def read_frame(self, i_frame):
        return self.read_frames([i_frame])[0]

    def read_frames(self, i_frames: [int]) -> [Frame]:
        frames = []
        for i in i_frames:
            # skip to frame i and read only this frame
            self._file.seek(self._frame_offset[i])
            string = self._file.read(self._frame_offset[i + 1] - self._frame_offset[i])
            frames.append(self.read_frame_from_string(string))

        return frames

    def read_frame_from_string(self, string: str):
        frame = Frame(self.n_atom)
        lines = string.splitlines()
        frame.step = int(lines[1])
        frame.xlo, frame.xhi = tuple(map(float, lines[5].split()))
        frame.ylo, frame.yhi = tuple(map(float, lines[6].split()))
        frame.zlo, frame.zhi = tuple(map(float, lines[7].split()))
        frame.box = np.array([frame.xhi - frame.xlo, frame.yhi - frame.ylo, frame.zhi - frame.zlo])
        tokens = lines[8].split()
        title = tokens[2:]

        data = '\n'.join(lines[9:])
        df = pd.read_csv(StringIO(data), header=None, index_col=None, names=title, sep='\s+')
        if 'ix' in df.columns:
            df.x += df['ix'] * frame.box[0]  # ix is revered words for pandas, so use df['ix'] instead of df.ix
            df.y += df['iy'] * frame.box[1]
            df.z += df['iz'] * frame.box[2]
        for row in df.itertuples():
            frame.charges[row.id-1] = float(row.q)
            if 'x' in df.columns:
                frame.positions[row.id-1] = np.array([row.x, row.y, row.z])
            elif 'xu' in df.columns:
                frame.positions[row.id-1] = np.array([row.xu, row.yu, row.zu])

        return frame

    def close(self):
        self._file.close()
