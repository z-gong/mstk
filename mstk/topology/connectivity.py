import numpy as np
from .atom import Atom
from .geometry import periodic_distance, periodic_angle, periodic_dihedral

__all__ = [
    'Bond',
    'Angle',
    'Dihedral',
    'Improper',
]


class Bond():
    '''
    A bond between two atoms.

    Unlike BondTerm in ForceField, the atoms are not sorted in connectivity
    because the properties of atoms (name, type, etc...) are prone to changing.
    It's better to sort them when compare bonds or match force field parameters.

    Parameters
    ----------
    atom1 : Atom
    atom2 : Atom
    order : int, Optional

    Attributes
    ----------
    atom1 : Atom
    atom2 : Atom
    order : int
        The integer bond order. See :class:`Bond.Order`
    '''

    class Order:
        '''
        Enumerator of integer bond orders

        The bonds in aromatic rings are in kekulized form.
        '''

        #: unspecified bond order
        UNSPECIFIED = 0
        #: order of single bond
        SINGLE = 1
        #: order of double bond
        DOUBLE = 2
        #: order of triple bond
        TRIPLE = 3

    def __init__(self, atom1, atom2, order=Order.UNSPECIFIED):
        self.atom1 = atom1
        self.atom2 = atom2
        self._order = order

    def __repr__(self):
        return f'<Bond: {self.name} {self.order}>'

    def equals(self, other):
        '''
        Check if two bonds represent the same connectivity.
        Return True if they contain the identical atoms, regardless of the sequence and bond order.

        Parameters
        ----------
        other : Bond

        Returns
        -------
        equal : bool
        '''
        if type(other) != Bond:
            return False
        return {self.atom1, self.atom2} == {other.atom1, other.atom2}

    @property
    def order(self):
        return self._order

    @order.setter
    def order(self, val):
        if val != self._order and self.atom1.molecule is not None:
            self.atom1.molecule._is_rdmol_valid = False
        self._order = val

    @property
    def name(self):
        '''
        Name of this bond

        Returns
        -------
        name : str
        '''
        return '%s-%s' % (self.atom1.name, self.atom2.name)

    @property
    def atoms(self):
        '''
        The atoms forming this bond

        Returns
        -------
        atoms : tuple of Atom
        '''
        return self.atom1, self.atom2

    @property
    def is_drude(self):
        '''
        Whether or not this bond is a Drude dipole bond

        Returns
        -------
        is : bool
        '''
        return self.atom1.is_drude or self.atom2.is_drude

    def evaluate(self, cell=None):
        '''
        Evaluate the length of this bond

        Parameters
        ----------
        box : UnitCell, Optional

        Returns
        -------
        value : float
        '''
        if cell:
            return periodic_distance(self.atom1.position, self.atom2.position, cell.get_size())

        delta = self.atom2.position - self.atom1.position
        return float(np.sqrt(delta.dot(delta)))

    @property
    def rdbond(self):
        '''
        The `rdkit.Chem.Bond` object associated with this bond.

        Returns
        -------
        rdbond : rdkit.Chem.Bond
        '''
        return self.atom1.molecule.rdmol.GetBondBetweenAtoms(self.atom1.id_in_mol, self.atom2.id_in_mol)

    @property
    def is_aromatic(self):
        return self.rdbond.GetIsAromatic()

    @property
    def is_conjugated(self):
        return self.rdbond.GetIsConjugated()

    @property
    def is_rotatable(self):
        if self.order != Bond.Order.SINGLE:
            return False
        if len(self.atom1.bonds) == 1 or (len(self.atom2.bonds)) == 1:
            return False

        return not self.rdbond.IsInRing()


class Angle():
    '''
    A angle between three atoms. The second atom is the central atom.

    Parameters
    ----------
    atom1 : Atom
    atom2 : Atom
    atom3 : Atom

    Attributes
    ----------
    atom1 : Atom
    atom2 : Atom
    atom3 : Atom
    '''

    def __init__(self, atom1, atom2, atom3):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3

    def __repr__(self):
        return '<Angle: %s>' % self.name

    def equals(self, other):
        '''
        Check if two angles represent the same connectivity.
        Return True if they contain the identical side atoms and center atom, regardless of the sequence of side atoms.

        Parameters
        ----------
        other : Angle

        Returns
        -------
        equal : bool
        '''
        if type(other) != Angle:
            return False
        if self.atom2 != other.atom2:
            return False
        return {self.atom1, self.atom3} == {other.atom1, other.atom3}

    @property
    def name(self):
        '''
        Name of this angle

        Returns
        -------
        name : str
        '''
        return '%s-%s-%s' % (self.atom1.name, self.atom2.name, self.atom3.name)

    @property
    def atoms(self):
        '''
        The atoms forming this angle

        Returns
        -------
        atoms : tuple of Atom
        '''
        return self.atom1, self.atom2, self.atom3

    @property
    def bonds(self):
        '''
        The two bonds between each side atom and the center atom

        Returns
        -------
        bonds : tuple of Bond
        '''
        _bond12 = Bond(self.atom1, self.atom2)
        _bond23 = Bond(self.atom2, self.atom3)
        try:
            bond12 = next(b for b in self.atom2.bonds if b.equals(_bond12))
        except StopIteration:
            raise Exception(f'{_bond12} not found in {self.atom2.molecule}')
        try:
            bond23 = next(b for b in self.atom2.bonds if b.equals(_bond23))
        except StopIteration:
            raise Exception(f'{_bond23} not found in {self.atom2.molecule}')

        return bond12, bond23

    def evaluate(self, cell=None):
        '''
        Evaluate the value of this angle in unit of radian

        Parameters
        ----------
        cell : UnitCell, Optional

        Returns
        -------
        value : float
        '''
        if cell:
            return periodic_angle(self.atom1.position, self.atom2.position, self.atom3.position, cell.get_size())

        vec1 = self.atom1.position - self.atom2.position
        vec2 = self.atom3.position - self.atom2.position
        cos = vec1.dot(vec2) / np.sqrt(vec1.dot(vec1) * vec2.dot(vec2))
        return float(np.arccos(np.clip(cos, -1, 1)))


class Dihedral():
    '''
    A dihedral between four atoms.

    Parameters
    ----------
    atom1 : Atom
    atom2 : Atom
    atom3 : Atom
    atom4 : Atom

    Attributes
    ----------
    atom1 : Atom
    atom2 : Atom
    atom3 : Atom
    atom4 : Atom
    override_atom_types: list of str
        When matching dihedral terms from FF, the type attributes of the four atoms will be overridden with these specified types.
        Each element can be a string or None. If None, the type attribute of the corresponding atom will be used, and not get overridden.
    '''

    def __init__(self, atom1, atom2, atom3, atom4):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.override_atom_types = [None, None, None, None]

    def __repr__(self):
        return '<Dihedral: %s>' % self.name

    def equals(self, other):
        '''
        Check if two dihedrals represent the same connectivity.
        Dihedrals with reversed sequence are considered as equal. E.g., i-j-k-l and l-k-j-i are the same.

        Parameters
        ----------
        other : Dihedral

        Returns
        -------
        equal : bool
        '''
        if type(other) != Dihedral:
            return False
        return (self.atom1 == other.atom1 and self.atom2 == other.atom2 and
                self.atom3 == other.atom3 and self.atom4 == other.atom4) \
               or (self.atom1 == other.atom4 and self.atom2 == other.atom3 and
                   self.atom3 == other.atom2 and self.atom4 == other.atom1)

    @property
    def name(self) -> str:
        '''
        Name of this dihedral

        Returns
        -------
        name : str
        '''
        return '%s-%s-%s-%s' \
               % (self.atom1.name, self.atom2.name, self.atom3.name, self.atom4.name)

    @property
    def atoms(self):
        '''
        The atoms forming this dihedral

        Returns
        -------
        atoms : tuple of Atom
        '''
        return self.atom1, self.atom2, self.atom3, self.atom4

    @property
    def bonds(self):
        '''
        The three bonds in this dihedral

        Returns
        -------
        bonds : tuple of Bond
        '''
        _bond12 = Bond(self.atom1, self.atom2)
        _bond23 = Bond(self.atom2, self.atom3)
        _bond34 = Bond(self.atom3, self.atom4)
        mol = self.atom1.molecule
        try:
            bond12 = next(b for b in self.atom1.bonds if b.equals(_bond12))
        except StopIteration:
            raise Exception(f'{_bond12} not found in {mol}')
        try:
            bond23 = next(b for b in self.atom2.bonds if b.equals(_bond23))
        except StopIteration:
            raise Exception(f'{_bond23} not found in {mol}')
        try:
            bond34 = next(b for b in self.atom3.bonds if b.equals(_bond34))
        except StopIteration:
            raise Exception(f'{_bond34} not found in {mol}')

        return bond12, bond23, bond34

    @property
    def angles(self):
        '''
        The two angles in this dihedral

        Returns
        -------
        angles : tuple of Angle
        '''
        _angle123 = Angle(self.atom1, self.atom2, self.atom3)
        _angle234 = Angle(self.atom2, self.atom3, self.atom4)
        mol = self.atom1.molecule
        try:
            angle123 = next(a for a in mol.angles if a.equals(_angle123))
        except StopIteration:
            raise Exception(f'{_angle123} not found in {mol}')
        try:
            angle234 = next(a for a in mol.angles if a.equals(_angle234))
        except StopIteration:
            raise Exception(f'{_angle234} not found in {mol}')

        return angle123, angle234

    def evaluate(self, cell=None):
        '''
        Evaluate the value of this dihedral in unit of radian

        Parameters
        ----------
        cell : UnitCell, Optional

        Returns
        -------
        value : float
        '''
        if cell:
            return periodic_dihedral(self.atom1.position, self.atom2.position, self.atom3.position, self.atom4.position,
                                     cell.get_size())

        vec1 = self.atom2.position - self.atom1.position
        vec2 = self.atom3.position - self.atom2.position
        vec3 = self.atom4.position - self.atom3.position
        n1 = np.cross(vec1, vec2)
        n2 = np.cross(vec2, vec3)
        cos = n1.dot(n2) / np.sqrt(n1.dot(n1) * n2.dot(n2))
        value = float(np.arccos(np.clip(cos, -1, 1)))
        sign = 1 if vec1.dot(n2) >= 0 else -1
        return sign * value


class Improper():
    '''
    An improper between four atoms. The first atom is the central atom.

    Parameters
    ----------
    atom1 : Atom
    atom2 : Atom
    atom3 : Atom
    atom4 : Atom

    Attributes
    ----------
    atom1 : Atom
    atom2 : Atom
    atom3 : Atom
    atom4 : Atom
    '''

    def __init__(self, atom1, atom2, atom3, atom4):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

    def __repr__(self):
        return '<Improper: %s>' % self.name

    def equals(self, other):
        '''
        Check if two impropers represent the same connectivity.
        Impropers are considered as equal if they have the same side atoms and central atom, regardless of the sequence of side atoms.
        So i-j-k-l and i-k-l-j are the same.

        Parameters
        ----------
        other : Improper

        Returns
        -------
        equal : bool
        '''
        if type(other) != Improper:
            return False
        if self.atom1 != other.atom1:
            return False
        return {self.atom2, self.atom3, self.atom4} == {other.atom2, other.atom3, other.atom4}

    @property
    def name(self):
        '''
        Name of this improper

        Returns
        -------
        name : str
        '''
        return '%s-%s-%s-%s' % (self.atom1.name, self.atom2.name, self.atom3.name, self.atom4.name)

    @property
    def atoms(self):
        '''
        The atoms forming this improper torsion

        Returns
        -------
        atoms : tuple of Atom
        '''
        return self.atom1, self.atom2, self.atom3, self.atom4

    def evaluate(self, cell=None):
        '''
        Evaluate the value of this improper torsion in unit of radian.
        The improper is defined as the angle between plane a1-a2-a3 and a2-a3-a4.

        Parameters
        ----------
        cell : UnitCell, Optional

        Returns
        -------
        value : float
        '''
        if cell:
            return periodic_dihedral(self.atom1.position, self.atom2.position, self.atom3.position, self.atom4.position,
                                     cell.get_size())

        vec1 = self.atom2.position - self.atom1.position
        vec2 = self.atom3.position - self.atom2.position
        vec3 = self.atom4.position - self.atom3.position
        n1 = np.cross(vec1, vec2)
        n2 = np.cross(vec2, vec3)
        cos = n1.dot(n2) / np.sqrt(n1.dot(n1) * n2.dot(n2))
        value = float(np.arccos(np.clip(cos, -1, 1)))
        sign = 1 if vec1.dot(n2) >= 0 else -1
        return sign * value
