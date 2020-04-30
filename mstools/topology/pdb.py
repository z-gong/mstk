import math
from ..topology import Topology
from ..forcefield import Element
from .atom import Atom
from .molecule import Molecule
from .topology import Topology


class Pdb(Topology):
    def __init__(self, file, **kwargs):
        super().__init__()
        self.parse(file)

    def parse(self, file):
        '''
        TODO To be implemented
        '''
        raise Exception('Reading PDB topology haven\'t been implemented')

    @staticmethod
    def save_to(top, file):
        if not top.has_position:
            raise Exception('Position is required for writing PDB file')

        lengths = top.cell.lengths * 10
        angles = top.cell.angles
        with open(file, 'w') as f:
            f.write('REMARK   Created by mstools\n')
            f.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 \n' % (
                lengths[0], lengths[1], lengths[2], angles[0], angles[1], angles[2]))

            for atom in top.atoms:
                pos = atom.position * 10  # convert from nm to A
                line = 'HETATM%5d %4s %4s %4d    %8.3f%8.3f%8.3f                      %2s\n' % (
                    (atom.id + 1) % 100000, atom.name[:4], atom.molecule.name[:4],
                    (atom.molecule.id + 1) % 10000, pos[0], pos[1], pos[2], atom.symbol[:2])
                f.write(line)

            for atom in top.atoms:
                partners = [a for a in atom.bond_partners if a.id > atom.id]
                if len(partners) > 0:
                    line = 'CONECT%5d' % (atom.id + 1)
                    for i, partner in enumerate(partners):
                        if i > 0 and i % 4 == 0:
                            line += '\nCONECT%5d' % (atom.id + 1)
                        line += '%5d' % (partner.id + 1)
                    line += '\n'
                    f.write(line)
