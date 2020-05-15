from .atom import Atom
from .molecule import Molecule
from .topology import Topology
from ..forcefield import Element


class XyzTopology(Topology):
    '''
    XYZ format only records the type and position of atoms
    There is no real topology, so all atoms are assumed to be in the same molecule
    The first column is treated as atom type instead of name or symbol
    '''

    def __init__(self, file, **kwargs):
        super().__init__()
        self.parse(file)

    def parse(self, file):
        with open(file) as f:
            n_atom = int(f.readline().strip())
            mol = Molecule()
            try:
                mol.name = f.readline().strip().split()[0]
            except:
                pass
            for i in range(n_atom):
                words = f.readline().split()
                atom = Atom()
                atom.type = words[0]
                element = Element.guess_from_atom_type(atom.type)
                atom.symbol = element.symbol
                atom.mass = element.mass
                atom.name = atom.symbol + str(mol.n_atom + 1)
                atom.position = tuple(map(lambda x: float(x) / 10, words[1:4]))
                mol.add_atom(atom)

        self.update_molecules([mol])

    @staticmethod
    def save_to(top, file):
        if not top.has_position:
            raise Exception('Position is required for writing XYZ file')

        with open(file, 'w')as f:
            f.write('%i\n' % top.n_atom)
            f.write('%s\n' % top.molecules[0].name)
            for atom in top.atoms:
                pos = atom.position * 10
                f.write('%-8s %11.5f %11.5f %11.5f\n' % (
                    atom.type or atom.symbol, pos[0], pos[1], pos[2]))
