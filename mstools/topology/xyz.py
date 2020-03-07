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
            raise Exception('Writing is not supported for XYZTopology')

    def parse(self):
        n_atom = int(self._file.readline().strip())
        self.remark = self._file.readline().strip()

        mol = Molecule()
        for i in range(n_atom):
            words = self._file.readline().split()
            atom = Atom()
            atom.type = words[0]
            mol.add_atom(atom)

        self.init_from_molecules([mol])

        self._file.close()
