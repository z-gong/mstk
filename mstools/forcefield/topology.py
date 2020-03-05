import numpy as np


class Atom():
    def __init__(self, name='UNK'):
        self.id = -1
        self.molecule: Molecule = None
        self.name = name
        self.type = ''
        self.symbol = ''
        self.charge = 0.
        self.mass = 0.
        self.position = np.array([0, 0, 0], dtype=float)

    def __repr__(self):
        return f'<Atom: {self.name} {self.id} {self.type}>'

    def __lt__(self, other):
        return self.name < other.name

    def __gt__(self, other):
        return self.name > other.name


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


class Bond():
    def __init__(self, atom1: Atom, atom2: Atom):
        self.atom1 = atom1
        self.atom2 = atom2

    def __eq__(self, other):
        return (self.atom1 == other.atom1 and self.atom2 == other.atom2) \
               or (self.atom1 == other.atom2 and self.atom2 == other.atom1)


class Angle():
    def __init__(self, atom1: Atom, atom2: Atom, atom3: Atom):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3

    def __eq__(self, other):
        if self.atom2 != other.atom2:
            return False

        return (self.atom1 == other.atom1 and self.atom3 == other.atom3) \
               or (self.atom1 == other.atom3 and self.atom3 == other.atom1)


class Dihedral():
    def __init__(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

    def __eq__(self, other):
        return (self.atom1 == other.atom1 and self.atom2 == other.atom2 and
                self.atom3 == other.atom3 and self.atom4 == other.atom4) \
               or (self.atom1 == other.atom4 and self.atom2 == other.atom3 and
                   self.atom3 == other.atom2 and self.atom4 == other.atom1)


class Improper():
    '''
    center atom is the first
    '''

    def __init__(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

    def __eq__(self, other):
        if self.atom1 != other.atom1:
            return False

        at12, at13, at14 = sorted([self.atom2, self.atom3, self.atom4])
        at22, at23, at24 = sorted([other.atom2, other.atom3, other.atom4])
        return at12 == at22 and at13 == at23 and at14 == at24


class Topology():
    def __init__(self):
        self.n_atom = 0
        self.atoms: [Atom] = []
        self.n_molecule = 0
        self.molecules: [Molecule] = []
        self.remark = ''
        self.is_drude = False

    def init_from_molecules(self, molecules: [Molecule]):
        self.n_molecule = len(molecules)
        self.molecules: [Molecule] = molecules[:]
        self.assign_id()

    def assign_id(self):
        idx_atom = 0
        for i, mol in enumerate(self.molecules):
            mol.id = i
            for j, atom in mol.atoms:
                atom.id = idx_atom
                idx_atom += 1

    @staticmethod
    def open(file, mode='r'):
        from .psf import PSF
        from .lammps import LammpsData
        from .xyz import XYZTopology

        if file.endswith('.psf'):
            return PSF(file, mode)
        elif file.endswith('.lmp'):
            return LammpsData(file, mode)
        elif file.endswith('.xyz'):
            return XYZTopology(file, mode)
        else:
            raise Exception('filename for topology not understand')
