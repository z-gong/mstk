import numpy as np


class Atom():
    def __init__(self, name='UNK'):
        self.id = -1
        self.name = name
        self.type = ''
        self.symbol = ''
        self.charge = 0.
        self.mass = 0.
        self.position = np.array([0, 0, 0], dtype=float)
        self._molecule: Molecule = None
        self._bonds: [Bond] = []

    def __repr__(self):
        return f'<Atom: {self.name} {self.id} {self.type}>'

    def __lt__(self, other):
        return self.name < other.name

    def __gt__(self, other):
        return self.name > other.name

    @property
    def molecule(self):
        return self._molecule

    @property
    def bonds(self):
        return self._bonds

    @property
    def bond_partners(self):
        return [bond.atom2 if bond.atom1 == self else bond.atom1 for bond in self._bonds]


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


class Molecule():
    def __init__(self, name='UNK'):
        self.id = -1
        self.name = name
        self._atoms = []
        self._bonds = []

    def add_atom(self, atom: Atom):
        self.atoms.append(atom)
        atom._molecule = self

    def remove_atom(self, atom: Atom):
        self.atoms.remove(atom)
        atom._molecule = None

    def add_bond(self, atom1: Atom, atom2: Atom):
        bond = Bond(atom1, atom2)
        self._bonds.append(bond)
        atom1._bonds.append(bond)
        atom2._bonds.append(bond)

    def remove_bond(self, bond: Bond):
        self._bonds.remove(bond)
        bond.atom1._bonds.remove(bond)
        bond.atom2._bonds.remove(bond)

    @property
    def atoms(self):
        return self._atoms

    @property
    def bonds(self):
        return self._bonds

    @property
    def initiated(self):
        return self.id != -1

    def __repr__(self):
        return f'<Molecule: {self.name} {self.id}>'


class Topology():
    def __init__(self):
        self.remark = ''
        self.is_drude = False
        self._atoms: [Atom] = []
        self._molecules: [Molecule] = []

    def init_from_molecules(self, molecules: [Molecule]):
        '''
        initialize a topology from a bunch of molecules
        note that there's no deepcopy
        the molecule and atoms are passed as reference
        '''
        self._molecules = molecules[:]
        self._atoms = [atom for mol in molecules for atom in mol.atoms]
        self.assign_id()

    def assign_id(self):
        idx_atom = 0
        for i, mol in enumerate(self._molecules):
            mol.id = i
            for j, atom in enumerate(mol.atoms):
                atom.id = idx_atom
                idx_atom += 1

    def add_molecule(self, molecule: Molecule):
        molecule.id = self.n_molecule
        self._molecules.append(molecule)
        for atom in molecule.atoms:
            atom.id = self.n_atom
            self._atoms.append(atom)

    @property
    def n_molecule(self):
        return len(self._molecules)

    @property
    def n_atom(self):
        return len(self._atoms)

    @property
    def molecules(self):
        return self._molecules

    @property
    def atoms(self):
        return self._atoms

    @staticmethod
    def open(file, mode='r'):
        from .psf import Psf
        from .lammps import LammpsData
        from .xyz import XyzTopology

        if file.endswith('.psf'):
            return Psf(file, mode)
        elif file.endswith('.lmp'):
            return LammpsData(file, mode)
        elif file.endswith('.xyz'):
            return XyzTopology(file, mode)
        else:
            raise Exception('filename for topology not understand')
