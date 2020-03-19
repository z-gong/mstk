import numpy as np
import pandas as pd
from io import StringIO
from . import Trajectory, Frame


class LammpsTrj(Trajectory):
    '''
    Read step, box and atomic positions (and charges optionally) from dump file of LAMMPS
    velocities will be ignored
    Because the topology information are detailed in data file, the mol, type, element in dump file will be ignored
    The length unit will be converted from A to nm
    '''

    def __init__(self, trj_file, mode='r'):
        super().__init__()
        self._file = open(trj_file, mode)
        if mode == 'r':
            self._get_info()
        elif mode == 'w':
            raise Exception('Writing support for LammpsTrj haven\'t been implemented')

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

        self._frame = Frame(self.n_atom)

    def read_frame(self, i_frame):
        # skip to frame i and read only this frame
        self._file.seek(self._frame_offset[i_frame])
        string = self._file.read(self._frame_offset[i_frame + 1] - self._frame_offset[i_frame])
        self._read_frame_from_string(string, self._frame)

        return self._frame

    def read_frames(self, i_frames: [int]) -> [Frame]:
        frames = []
        for i in i_frames:
            frame = Frame(self.n_atom)
            # skip to frame i and read only this frame
            self._file.seek(self._frame_offset[i])
            string = self._file.read(self._frame_offset[i + 1] - self._frame_offset[i])
            self._read_frame_from_string(string, frame)
            frames.append(frame)

        return frames

    def _read_frame_from_string(self, string: str, frame: Frame):
        lines = string.splitlines()
        frame.step = int(lines[1])
        try:
            xlo, xhi = tuple(map(lambda x: float(x) / 10, lines[5].split()))  # convert from A to nm
            ylo, yhi = tuple(map(lambda x: float(x) / 10, lines[6].split()))
            zlo, zhi = tuple(map(lambda x: float(x) / 10, lines[7].split()))
            frame.cell.set_box([xhi - xlo, yhi - ylo, zhi - zlo])
        except:
            # read triclinic box, convert from A to nm
            xlo, xhi, bx = tuple(map(lambda x: float(x) / 10, lines[5].split()))
            ylo, yhi, cx = tuple(map(lambda x: float(x) / 10, lines[6].split()))
            zlo, zhi, cy = tuple(map(lambda x: float(x) / 10, lines[7].split()))
            frame.cell.set_box([[xhi - xlo, 0, 0],[bx,  yhi - ylo, 0],[cx, cy ,zhi - zlo]])

        tokens = lines[8].split()
        title = tokens[2:]

        data = '\n'.join(lines[9:])
        df = pd.read_csv(StringIO(data), header=None, index_col=None, names=title, sep=r'\s+')
        wrapped = False
        if 'x' in df.columns:
            df['x'] = df['x'] / 10 - xlo  # convert from A to nm
            df['y'] = df['y'] / 10 - ylo
            df['z'] = df['z'] / 10 - zlo
            wrapped = True
        elif 'xs' in df.columns:
            df['x'] = df.xs * frame.cell.size[0]
            df['y'] = df.ys * frame.cell.size[1]
            df['z'] = df.zs * frame.cell.size[2]
            wrapped = True
        elif 'xu' in df.columns:
            df['x'] = df.xu / 10 - xlo  # convert from A to nm
            df['y'] = df.yu / 10 - ylo
            df['z'] = df.zu / 10 - zlo
        elif 'xsu' in df.columns:
            df['x'] = df.xsu * frame.cell.size[0]
            df['y'] = df.ysu * frame.cell.size[1]
            df['z'] = df.zsu * frame.cell.size[2]
        if wrapped:
            if 'ix' not in df.columns:
                print('warning: image flag not found for wrapped positions')
            else:
                # ix is revered words for pandas, so use df['ix'] instead of df.ix
                df.x += df['ix'] * frame.cell.size[0]
                df.y += df['iy'] * frame.cell.size[1]
                df.z += df['iz'] * frame.cell.size[2]

        frame.has_charge = 'q' in df.columns
        for row in df.itertuples():
            frame.positions[row.id - 1] = np.array([row.x, row.y, row.z])
            if frame.has_charge:
                frame.charges[row.id - 1] = float(row.q)

        return frame
