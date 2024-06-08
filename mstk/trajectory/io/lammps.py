import numpy as np
import pandas as pd
from io import StringIO
from mstk import logger
from mstk.trajectory.handler import TrjHandler


class LammpsTrj(TrjHandler):
    '''
    Read step, cell and atomic positions (and charges if provided) from dump file of Lammps.

    Velocities are ignored.
    Because the topology information are detailed in data file, the mol, type, element in dump file are also ignored.
    The cell used by Lammps is not originated at (0, 0, 0) usually.
    This handler will translate the positions of atoms based on the lower bound of the cell.

    Writing to LammpsDump is not supported yet.
    '''

    def __init__(self, trj_file, mode='r'):
        super().__init__()

        if mode == 'r':
            # open it in binary mode so that we can correctly seek despite of line ending
            self._file = open(trj_file, 'rb')
        else:
            raise Exception('Writing support for LammpsTrj haven\'t been implemented')

    def get_info(self):
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
        # the last element is the length of whole file
        self._line_offset = [0]
        offset = 0
        for line in self._file:
            offset += len(line)
            self._line_offset.append(offset)
        self._file.seek(0)

        # build a list of frame offsets
        self.n_frame = len(self._line_offset) // (9 + self.n_atom)
        self._frame_offset = []
        for i in range(self.n_frame + 1):
            line_start = (9 + self.n_atom) * i
            self._frame_offset.append(self._line_offset[line_start])

        return self.n_atom, self.n_frame

    def read_frame(self, i_frame, frame):
        # skip to frame i and read only this frame
        self._file.seek(self._frame_offset[i_frame])
        string = self._file.read(self._frame_offset[i_frame + 1] - self._frame_offset[i_frame]).decode()
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
            frame.cell.set_box([[xhi - xlo, 0, 0], [bx, yhi - ylo, 0], [cx, cy, zhi - zlo]])
        box = frame.cell.get_size()

        tokens = lines[8].split()
        title = tokens[2:]

        data = '\n'.join(lines[9:])
        df = pd.read_csv(StringIO(data), header=None, index_col=None, names=title, sep=r'\s+')
        wrapped = False
        if 'x' in df.columns:
            df['x'] = df['x'] / 10  # convert from A to nm
            df['y'] = df['y'] / 10
            df['z'] = df['z'] / 10
            wrapped = True
        elif 'xs' in df.columns:
            df['x'] = df.xs * box[0] + xlo
            df['y'] = df.ys * box[1] + ylo
            df['z'] = df.zs * box[2] + zlo
            wrapped = True
        elif 'xu' in df.columns:
            df['x'] = df.xu / 10  # convert from A to nm
            df['y'] = df.yu / 10
            df['z'] = df.zu / 10
        elif 'xsu' in df.columns:
            df['x'] = df.xsu * box[0] + xlo
            df['y'] = df.ysu * box[1] + ylo
            df['z'] = df.zsu * box[2] + zlo
        if wrapped:
            if 'ix' not in df.columns:
                logger.warning('Image flag not found for wrapped positions')
            else:
                # ix is revered words for pandas, so use df['ix'] instead of df.ix
                df['x'] += df['ix'] * box[0]
                df['y'] += df['iy'] * box[1]
                df['z'] += df['iz'] * box[2]

        frame.has_charge = 'q' in df.columns
        for row in df.itertuples():
            frame.positions[row.id - 1] = np.array([row.x, row.y, row.z])
            if frame.has_charge:
                frame.charges[row.id - 1] = float(row.q)


TrjHandler.register_format('.lammpstrj', LammpsTrj)
TrjHandler.register_format('.ltrj', LammpsTrj)
