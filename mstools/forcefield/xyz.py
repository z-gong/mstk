from . import Topology, Atom, Molecule

class XYZTopology(Topology):
    '''
    xyz format only records the type and position of atoms
    there is no real topology
    we treat the first column as type instead of name
    '''

    def __init__(self, file, mode='r'):
        super(XYZTopology, self).__init__()
        self._file = open(file, mode)
        if mode == 'r':
            self.parse()
        elif mode == 'w':
            raise Exception('Writing is not supported for XYZTopology')

    def parse(self):
        self.n_atom = int(self._file.readline().strip())
        self.remark = self._file.readline().strip()

        mol = Molecule()
        mol.id = 0
        self.molecules.append(mol)
        for i in range(self.n_atom):
            words = self._file.readline().split()
            atom = Atom()
            atom.id = i
            atom.type = words[0]
            mol.add_atom(atom)
            self.atoms.append(atom)

        self._file.close()
