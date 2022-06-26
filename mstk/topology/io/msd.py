import io
from mstk.topology.atom import Atom
from mstk.topology.molecule import Molecule
from mstk.topology.topology import Topology
from mstk.chem.element import Element


class Msd():
    '''
    Generate Topology from MSD file of DFF.

    The atom names, symbols, types, charges, positions, molecules, bonds and unit cell (if provided in the file) are parsed.
    The bond orders, formal charges and subsets are ignored.

    Parameters
    ----------
    file : str or file-like object
        MSD file
    kwargs : dict
        Ignored

    Attributes
    ----------
    topology : Topology

    Examples
    --------
    >>> msd = Msd('input.msd')
    >>> topology = msd.topology

    '''

    def __init__(self, file, **kwargs):
        if type(file) is str:
            with open(file) as f:
                content = f.read()
        elif isinstance(file, io.IOBase):
            content = file.read()
        else:
            raise Exception('A filename or a file object required')

        self.topology = Topology()
        self._parse(content)

    def _parse(self, content):
        lines = content.splitlines()
        # Ignore the first three lines
        # The forth line could be cell information or atom number
        line = lines[3]
        if line.startswith('PBC'):
            words = line.strip().split()
            lengths = list(map(lambda x: float(x) / 10, words[1:4]))
            angles = list(map(float, words[4:7]))
            self.topology.cell.set_box([lengths, angles])
            n_header = 4
        else:
            n_header = 3

        n_atom = int(lines[n_header])
        atoms = [Atom() for i in range(n_atom)]
        mol_names = {}
        for i in range(n_atom):
            line = lines[n_header + 1 + i]
            words = line.strip().split()
            atom_id = int(words[0])  # index in MSD starts from 1
            atom = atoms[atom_id - 1]
            atom.id = atom_id - 1  # atom.id stars from 0
            atom.name = words[1]
            element = Element(int(words[2]))
            atom.symbol = element.symbol
            atom.mass = element.mass
            atom.type = words[3]
            atom.charge, x, y, z = list(map(float, words[5:9]))
            atom.position = [x / 10, y / 10, z / 10]

            mol_id = int(words[9])  # starts from 1
            mol_name = words[10]
            if mol_id not in mol_names:
                mol_names[mol_id] = mol_name
            atom._mol_id = mol_id  # temporary attribute for identifying molecules

        molecules = [Molecule(name) for id, name in sorted(mol_names.items())]
        for atom in atoms:
            mol = molecules[atom._mol_id - 1]
            mol.add_atom(atom, update_topology=False)

        n_bond = int(lines[n_header + n_atom + 1])
        for i in range(n_bond):
            line = lines[n_header + n_atom + 2 + i]
            words = line.strip().split()
            atom1 = atoms[int(words[0]) - 1]
            atom2 = atoms[int(words[1]) - 1]
            atom1.molecule.add_bond(atom1, atom2)

        self.topology.update_molecules(molecules)
        self.topology.generate_angle_dihedral_improper()


Topology.register_format('.msd', Msd)
