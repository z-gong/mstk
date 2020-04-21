from ..topology import Topology
from ..forcefield import Element
from .atom import Atom
from .molecule import Molecule
from .topology import Topology


class Pdb(Topology):
    '''
    The reading haven't been implemented yet
    Anyway, i don't like pdb format
    '''

    def __init__(self, file):
        super().__init__()
        self.parse(file)

    def parse(self, file):
        raise Exception('Reading PDB topology haven\'t been implemented')

    @staticmethod
    def save_to(top, file):
        if not top.has_position:
            raise Exception('Position is required for writing PDB file')

        lengths = top.cell.lengths * 10
        angles = top.cell.angles
        with open(file, 'w') as f:
            f.write('TITLE    PDB written by mstools\n')
            f.write('REMARK   %s\n' % top.remark.replace('\n', ''))
            f.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 \n' % (
                lengths[0], lengths[1], lengths[2], angles[0], angles[1], angles[2]))

            for id in range(top.n_atom):
                atom = top.atoms[id]
                pos = top.positions[id] * 10  # convert from nm to A
                line = 'ATOM  %5d %4s %4s%5d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n' % (
                    (id + 1) % 100000, atom.name[:4], atom.molecule.name[:4],
                    (atom.molecule.id + 1) % 100000, pos[0], pos[1], pos[2], atom.symbol)
                f.write(line)
