from .atom import Atom
from .molecule import Molecule
from .topology import Topology
from ..forcefield import Element


class XyzTopology(Topology):
    '''
    xyz format only records the type and position of atoms
    there is no real topology
    we treat the first column as type instead of name
    '''

    def __init__(self, file):
        super(XyzTopology, self).__init__()
        self.parse(file)

    def parse(self, file):
        with open(file) as f:
            n_atom = int(f.readline().strip())
            self.remark = f.readline().strip()

            mol = Molecule()
            for i in range(n_atom):
                words = f.readline().split()
                atom = Atom()
                atom.type = words[0]
                atom.symbol = Element.guess_from_atom_type(atom.type).symbol
                atom.name = atom.symbol + str(mol.n_atom + 1)
                atom.position = tuple(map(lambda x: float(x) / 10, words[1:4]))
                atom.has_position = True
                mol.add_atom(atom)

        self.init_from_molecules([mol])

    @staticmethod
    def save_to(top, file):
        if not top.has_position:
            raise Exception('Position is required for writing xyz file')

        with open(file, 'w')as f:
            f.write('%i\n' % top.n_atom)
            f.write('%s\n' % top.remark)
            for atom in top.atoms:
                pos = atom.position * 10
                f.write('%-8s %11.5f %11.5f %11.5f\n' % (atom.type, pos[0], pos[1], pos[2]))
