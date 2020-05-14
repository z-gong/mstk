from .atom import Atom


class Bond():
    def __init__(self, atom1: Atom, atom2: Atom):
        '''
        Because the properties of atoms are prone to changing,
        there's no point to sort atoms here
        it's better to sort them when compare bonds or match force field parameters
        '''
        self.atom1 = atom1
        self.atom2 = atom2

    def __repr__(self):
        return '<Bond: %s>' % self.name

    def __eq__(self, other):
        if type(other) != Bond:
            return False
        return {self.atom1, self.atom2} == {other.atom1, other.atom2}

    @property
    def name(self) -> str:
        return '%s-%s' % (self.atom1.name, self.atom2.name)

    @property
    def is_drude(self) -> bool:
        return self.atom1.is_drude or self.atom2.is_drude


class Angle():
    def __init__(self, atom1: Atom, atom2: Atom, atom3: Atom):
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
        return '%s-%s-%s' % (self.atom1.name, self.atom2.name, self.atom3.name)


class Dihedral():
    def __init__(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom):
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
        return '%s-%s-%s-%s' \
               % (self.atom1.name, self.atom2.name, self.atom3.name, self.atom4.name)


class Improper():
    '''
    center atom is the first
    '''

    def __init__(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom):
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
    def name(self) -> str:
        return '%s-%s-%s-%s' % (self.atom1.name, self.atom2.name, self.atom3.name, self.atom4.name)
