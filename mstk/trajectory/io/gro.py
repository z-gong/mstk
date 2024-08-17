import numpy as np
from mstk.topology import Topology
from mstk.trajectory import Frame
from mstk.trajectory.handler import TrjHandler


class Gro(TrjHandler):
    '''
    Read and write cell, atomic positions and optionally velocities from/to GRO file.
    '''

    def __init__(self, file, mode='r'):
        super().__init__()
        if mode not in ('r', 'w', 'a'):
            raise Exception('Invalid mode')

        if mode == 'r':
            # open it in binary mode so that we can correctly seek despite of line ending
            self._file = open(file, 'rb')
        elif mode == 'a':
            self._file = open(file, 'ab')
        elif mode == 'w':
            self._file = open(file, 'wb')

    def get_info(self):
        try:
            self._file.readline()
            self.n_atom = int(self._file.readline())
        except:
            print('Invalid gro file')
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
        self.n_frame = len(self._line_offset) // (3 + self.n_atom)
        self._frame_offset = []
        for i in range(self.n_frame + 1):
            line_start = (3 + self.n_atom) * i
            self._frame_offset.append(self._line_offset[line_start])

        return self.n_atom, self.n_frame

    def read_frame(self, i_frame, frame):
        # skip to frame i and read only this frame
        self._file.seek(self._frame_offset[i_frame])
        lines = self._file.read(self._frame_offset[i_frame + 1] - self._frame_offset[i_frame]) \
            .decode().splitlines()
        # assume there are velocities. we'll see later
        frame.has_velocity = True
        for i in range(self.n_atom):
            line = lines[i + 2]
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            frame.positions[i][:] = x, y, z

            if frame.has_velocity:
                try:
                    vx = float(line[44:52])
                    vy = float(line[52:60])
                    vz = float(line[60:68])
                except:
                    frame.has_velocity = False
                else:
                    frame.velocities[i][:] = vx, vy, vz
        _box = tuple(map(float, lines[self.n_atom + 2].split()))
        if len(_box) == 3:
            frame.cell.set_box(_box)
        elif len(_box) == 9:
            ax, by, cz, ay, az, bx, bz, cx, cy = _box
            frame.cell.set_box([[ax, ay, az], [bx, by, bz], [cx, cy, cz]])
        else:
            raise ValueError('Invalid box')

    def write_frame(self, frame, topology, subset=None, write_velocity=False, **kwargs):
        '''
        Write a frame into the opened GRO file

        Parameters
        ----------
        frame : Frame
        topology : Topology
        subset : list of int, optional
        write_velocity : bool
            Whether or not velocities should be written.
            If set to True but velocities not available in frame, an Exception will be raised.
        kwargs : dict
            Ignored
        '''
        if subset is None:
            subset = list(range(len(frame.positions)))
        if write_velocity and not frame.has_velocity:
            raise Exception('Velocities are requested but not exist in frame')

        string = 'Created by mstk: step= %i, t= %f ps\n' % (frame.step, frame.time)
        string += '%i\n' % len(subset)

        for id in subset:
            atom = topology.atoms[id]
            residue = atom.residue
            pos = frame.positions[id]
            if any(np.abs(pos) >= 1000):
                raise Exception('Positions are too large to be written in GRO format')
            string += '%5i%5s%5s%5i%8.3f%8.3f%8.3f' % (
                (residue.id + 1) % 100000, residue.name[:5], atom.symbol[:5], (atom.id + 1) % 100000,
                pos[0], pos[1], pos[2])
            if write_velocity:
                vel = frame.velocities[id]
                string += '%8.4f%8.4f%8.4f' % (vel[0], vel[1], vel[2])
            string += '\n'

        a, b, c = frame.cell.vectors
        string += ' %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n' % (
            a[0], b[1], c[2], a[1], a[2], b[0], b[2], c[0], c[1])

        self._file.write(string.encode())
        self._file.flush()

TrjHandler.register_format('.gro', Gro)
