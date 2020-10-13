from .atom import Atom


class Bond():
    '''
    A bond between two atoms.

    Two bonds are considered as equal if they contain the identical atoms, regardless of the sequence.

    Unlike BondTerm in ForceField, the atoms are not sorted in connectivity
    because the properties of atoms (name, type, etc...) are prone to changing.
    It's better to sort them when compare bonds or match force field parameters.

    Parameters
    ----------
    atom1 : Atom
    atom2 : Atom

    Attributes
    ----------
    atom1 : Atom
    atom2 : Atom
    '''
    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2

    def __repr__(self):
        return '<Bond: %s>' % self.name

    def __eq__(self, other):
        if type(other) != Bond:
            return False
        return {self.atom1, self.atom2} == {other.atom1, other.atom2}

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
    def is_drude(self):
        '''
        Whether or not this bond is a Drude dipole bond

        Returns
        -------
        is : bool
        '''
        return self.atom1.is_drude or self.atom2.is_drude


class Angle():
    '''
    A angle between three atoms.

    The second atom is the central atom.
    Two angles are considered as equal if they contain the identical side atoms and center atom,
    regardless of the sequence of side atoms.

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

    def __eq__(self, other):
        if type(other) != Angle:
            return False
        if self.atom2 != other.atom2:
            return False
        return {self.atom1, self.atom3} == {other.atom1, other.atom3}

    @property
    def name(self) -> str:
        '''
        Name of this angle

        Returns
        -------
        name : str
        '''
        return '%s-%s-%s' % (self.atom1.name, self.atom2.name, self.atom3.name)


class Dihedral():
    '''
    A dihedral between four atoms.

    Two dihedrals with reversed sequence are considered as equal.
    So i-j-k-l and l-k-j-i are the same.

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
        return '<Dihedral: %s>' % self.name

    def __eq__(self, other):
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


class Improper():
    '''
    A improper between four atoms.

    The first atom is the central atom.
    There's arbitrariness in the definition of improper.
    Herein, two impropers are considered as equal if they have the same side atoms and central atom,
    regardless of the sequence of side atoms.
    So i-j-k-l and i-k-l-j are the same.

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

    def __eq__(self, other):
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
