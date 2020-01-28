import traceback
import numpy as np
from io import StringIO
import pandas as pd
from . import Atom, Topology, Trajectory, Frame


class PDB(Trajectory):
    def __init__(self, file, mode='r'):
        super().__init__()
        self._file = open(file, mode)
        self._n_model = 0

    def write_frame(self, topology: Topology, frame: Frame, subset=None):
        self._file.write('TITLE    simulation box\n')
        self._file.write('REMARK   created by mstools\n')
        self._file.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 \n' % (
            frame.box[0], frame.box[1], frame.box[2], 90, 90, 90))
        self._n_model += 1
        self._file.write('MODEL    %i\n' % self._n_model)

        if subset is None:
            subset = list(range(topology.n_atom))
        for ii, id in enumerate(subset):
            atom = topology.atoms[id]
            position = frame.positions[id]
            line = 'ATOM  %5d %4s %4s %4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n' %(
                ii + 1, atom.element, atom.molecule.name[:4], atom.molecule.id + 1,
                position[0], position[1], position[2], atom.element)
            self._file.write(line)
        self._file.write('ENDMDL\n')

    def close(self):
        self._file.close()
