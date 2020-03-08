from ..forcefield import Element
from ..topology import Topology, Atom, Molecule


class XyzTopology(Topology):
    '''
    xyz format only records the type and position of atoms
    there is no real topology
    we treat the first column as type instead of name
    '''

    def __init__(self, file, mode='r'):
        super(XyzTopology, self).__init__()
        self._file = open(file, mode)
        if mode == 'r':
            self.parse()
        elif mode == 'w':
            pass

    def parse(self):
        n_atom = int(self._file.readline().strip())
        self.remark = self._file.readline().strip()

        mol = Molecule()
        for i in range(n_atom):
            words = self._file.readline().split()
            atom = Atom()
            atom.type = words[0]
            atom.symbol = Element.guess_from_atom_type(atom.type).symbol
            atom.name = atom.symbol + str(mol.n_atom + 1)
            atom.position = tuple(map(lambda x: float(x) / 10, words[1:4]))
            atom.has_position = True
            mol.add_atom(atom)

        self.init_from_molecules([mol])

    def write(self):
        if not self.has_position:
            raise Exception('Position is required for writing xyz file')
        self._file.write('%i\n' % self.n_atom)
        self._file.write('%s\n' % self.remark)
        for atom in self.atoms:
            pos = atom.position * 10
            self._file.write('%-8s %10.4f %10.4f %10.4f\n' % (
                atom.type, pos[0], pos[1], pos[2]))
