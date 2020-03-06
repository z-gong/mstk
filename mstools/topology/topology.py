import numpy as np


class Atom():
    def __init__(self, name='UNK'):
        self.id = -1
        self.name = name
        self.type = ''
        self.symbol = ''
        self.charge = 0.
        self.mass = 0.
        self.alpha = 0.  # for polarizable model
        self.thole = 0.  # for polarizable model
        self.has_position = False
        self.position = np.array([0, 0, 0], dtype=float)
        self.has_velocity = False
        self.velocity = np.array([0, 0, 0], dtype=float)
        self._molecule: Molecule = None
        self._bonds: [Bond] = []

    def __repr__(self):
        return f'<Atom: {self.name} {self.id} {self.type}>'

    def __lt__(self, other):
        return self.id < other.id

    def __gt__(self, other):
        return self.id > other.id

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
        self._angles = []
        self._dihedrals = []
        self._impropers = []

    def add_atom(self, atom: Atom):
        self.atoms.append(atom)
        atom._molecule = self

    def remove_atom(self, atom: Atom):
        self.atoms.remove(atom)
        atom._molecule = None

    def add_bond(self, atom1: Atom, atom2: Atom):
        # if atom1 not in self._atoms or atom2 not in self._atoms:
        #     raise Exception('both atoms %s %s should belong to same molecule %s'%(atom1, atom2, self))

        bond = Bond(atom1, atom2)
        self._bonds.append(bond)
        atom1._bonds.append(bond)
        atom2._bonds.append(bond)
        return bond

    def remove_bond(self, bond: Bond):
        self._bonds.remove(bond)
        bond.atom1._bonds.remove(bond)
        bond.atom2._bonds.remove(bond)

    def add_angle(self, atom1: Atom, atom2: Atom, atom3: Atom):
        # if atom1 not in self._atoms or atom2 not in self._atoms or atom3 not in self._atoms:
        #     raise Exception('both atoms %s %s %s should belong to same molecule %s' % (atom1, atom2, atom3, self))

        angle = Angle(atom1, atom2, atom3)
        self._angles.append(angle)
        return angle

    def remove_angle(self, angle: Angle):
        self._angles.remove(angle)

    def add_dihedral(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom):
        # if atom1 not in self._atoms or atom2 not in self._atoms \
        #     or atom3 not in self._atoms or atom4 not in self._atoms:
        #     raise Exception('both atoms %s %s %s %s should belong to same molecule %s' % (atom1, atom2, atom3, atom4, self))

        dihedral = Dihedral(atom1, atom2, atom3, atom4)
        self._dihedrals.append(dihedral)
        return dihedral

    def remove_dihedral(self, dihedral: Dihedral):
        self._dihedrals.remove(dihedral)

    def add_improper(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom):
        # if atom1 not in self._atoms or atom2 not in self._atoms \
        #         or atom3 not in self._atoms or atom4 not in self._atoms:
        #     raise Exception('both atoms %s %s %s %s should belong to same molecule %s' % (atom1, atom2, atom3, atom4, self))

        improper = Improper(atom1, atom2, atom3, atom4)
        self._impropers.append(improper)
        return improper

    def remove_improper(self, improper: Improper):
        self._impropers.remove(improper)

    @property
    def n_atom(self):
        return len(self._atoms)

    @property
    def n_bond(self):
        return len(self._bonds)

    @property
    def n_angle(self):
        return len(self._angles)

    @property
    def n_dihedral(self):
        return len(self._dihedrals)

    @property
    def n_improper(self):
        return len(self._impropers)

    @property
    def atoms(self):
        return self._atoms

    @property
    def bonds(self):
        return self._bonds

    @property
    def angles(self):
        return self._angles

    @property
    def dihedrals(self):
        return self._dihedrals

    @property
    def impropers(self):
        return self._impropers

    @property
    def initiated(self):
        return self.id != -1

    def __repr__(self):
        return f'<Molecule: {self.name} {self.id}>'

    def generate_angle_dihedral_improper(self):
        pass


class Topology():
    def __init__(self):
        self.remark = ''
        self.is_drude = False
        self._molecules: [Molecule] = []
        self._box = np.array([0, 0, 0], dtype=float)

    def init_from_molecules(self, molecules: [Molecule]):
        '''
        initialize a topology from a bunch of molecules
        note that there's no deepcopy
        the molecule and atoms are passed as reference
        '''
        self._molecules = molecules[:]
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

    def set_positions(self, positions):
        if self.n_atom != len(positions):
            raise Exception('Length of positions should equal to the number of atoms')
        for i, atom in enumerate(self.atoms):
            atom.position = np.array(positions[i])
            atom.has_position = True

    def set_velocities(self, velocities):
        if self.n_atom != len(velocities):
            raise Exception('Length of velocities should equal to the number of atoms')
        for i, atom in enumerate(self.atoms):
            atom.velocity = np.array(velocities[i])
            atom.has_velocity = True

    @property
    def n_molecule(self):
        return len(self._molecules)

    @property
    def n_atom(self):
        return sum([mol.n_atom for mol in self._molecules])

    @property
    def n_bond(self):
        return sum([mol.n_bond for mol in self._molecules])

    @property
    def n_angle(self):
        return sum([mol.n_angle for mol in self._molecules])

    @property
    def n_dihedral(self):
        return sum([mol.n_dihedral for mol in self._molecules])

    @property
    def n_improper(self):
        return sum([mol.n_improper for mol in self._molecules])

    @property
    def molecules(self):
        return self._molecules

    @property
    def atoms(self):
        return [atom for mol in self._molecules for atom in mol.atoms]

    @property
    def bonds(self):
        return [bond for mol in self._molecules for bond in mol.bonds]

    @property
    def angles(self):
        return [angle for mol in self._molecules for angle in mol.angles]

    @property
    def dihedrals(self):
        return [dihedral for mol in self._molecules for dihedral in mol.dihedrals]

    @property
    def impropers(self):
        return [improper for mol in self._molecules for improper in mol.impropers]

    @property
    def has_position(self):
        return all(atom.has_position for atom in self.atoms)

    @property
    def has_velocity(self):
        return all(atom.has_velocity for atom in self.atoms)

    @property
    def positions(self):
        return [atom.position for atom in self.atoms]

    @property
    def velocities(self):
        return [atom.velocity for atom in self.atoms]

    @property
    def box(self):
        return self._box

    @box.setter
    def box(self, value):
        if len(value) != 3:
            raise Exception('Invalid box size')
        self._box = np.array(value)

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
