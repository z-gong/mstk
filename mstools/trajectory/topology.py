import numpy as np


class Atom():
    def __init__(self, name='UNK'):
        self.id = -1
        self.molecule: Molecule = None
        self.name = name
        self.type = ''
        self.element = ''
        self.charge = 0.
        self.mass = 0.
        self.position = np.array([0, 0, 0], dtype=float)

    def __repr__(self):
        return f'<Atom: {self.name} {self.id} {self.type}>'


class Molecule():
    def __init__(self, name='UNK'):
        self.id = -1
        self.name = name
        self.atoms: [Atom] = []

    def add_atom(self, atom: Atom):
        self.atoms.append(atom)
        atom.molecule = self

    def remove_atom(self, atom: Atom):
        self.atoms.remove(atom)
        atom.molecule = None

    def __repr__(self):
        return f'<Molecule: {self.name} {self.id}>'

class Topology():
    def __init__(self):
        pass

    def init_from_molecules(self, molecules: [Molecule]):
        self.n_molecule = len(molecules)
        self.molecules:[Molecule] = molecules[:]
        self.assign_id()

    def assign_id(self):
        idx_atom = 0
        for i, mol in enumerate(self.molecules):
            mol.id = i
            for j, atom in mol.atoms:
                atom.id = idx_atom
                idx_atom += 1

atomic_nr = {'H': 1, 'Li': 3, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
             'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16,
             'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Ti': 22, 'Fe': 26,
             'Zn': 30, 'Se': 34, 'Br': 35, 'Kr': 36, 'Mo': 42, 'Ru': 44,
             'Sn': 50, 'Te': 52, 'I': 53, 'Xe': 54}

def atomic_symbol(name):
    if name[:2] in atomic_nr:
        return name[:2]
    elif name[0] in atomic_nr:
        return name[0]
    else:
        print('warning: unknown symbol for atom ' + name)
        return name

def atomic_number(name):
    return atomic_nr.get(atomic_symbol(name), 0)
