from ..topology import Topology
from . import Trajectory, Frame


class Pdb(Trajectory):
    '''
    Write frames into pdb file
    The reading haven't been implemented yet
    Anyway, i don't like pdb format
    '''

    def __init__(self, file, mode='r'):
        super().__init__()
        self._file = open(file, mode)
        if mode == 'r':
            raise Exception('Reading support for PDB haven\'t been implemented')
        elif mode == 'w':
            self._n_model = 0

    def write_frame(self, topology: Topology, frame: Frame, subset=None):
        self._file.write('TITLE    simulation box\n')
        self._file.write('REMARK   created by mstools\n')
        self._file.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 \n' % (
            frame.box[0] * 10, frame.box[1] * 10, frame.box[2] * 10, 90, 90, 90))
        self._n_model += 1
        self._file.write('MODEL    %i\n' % self._n_model)

        if subset is None:
            subset = list(range(topology.n_atom))
        for id in subset:
            atom = topology.atoms[id]
            pos = frame.positions[id] * 10  # convert from nm to A
            line = 'ATOM  %5d %4s %4s%5d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n' % (
                (id + 1) % 100000, atom.symbol, atom.molecule.name[:4], (atom.molecule.id + 1) % 100000,
                pos[0], pos[1], pos[2], atom.symbol)
            self._file.write(line)
        self._file.write('ENDMDL\n')

    def close(self):
        self._file.close()
