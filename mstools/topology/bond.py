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

    def __str__(self):
        return '<Bond: %s-%s>' % (self.atom1.name, self.atom2.name)

    def __repr__(self):
        return str(self) + ' at 0x' + str(hex(id(self))[2:].upper())

    def __eq__(self, other):
        return (self.atom1 == other.atom1 and self.atom2 == other.atom2) \
               or (self.atom1 == other.atom2 and self.atom2 == other.atom1)

    def __lt__(self, other):
        return [self.atom1, self.atom2] < [other.atom1, other.atom2]

    def __gt__(self, other):
        return [self.atom1, self.atom2] > [other.atom1, other.atom2]


class Angle():
    def __init__(self, atom1: Atom, atom2: Atom, atom3: Atom):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3

    def __str__(self):
        return '<Angle: %s-%s-%s>' % (self.atom1.name, self.atom2.name, self.atom3.name)

    def __repr__(self):
        return str(self) + ' at 0x' + str(hex(id(self))[2:].upper())

    def __eq__(self, other):
        if self.atom2 != other.atom2:
            return False

        return (self.atom1 == other.atom1 and self.atom3 == other.atom3) \
               or (self.atom1 == other.atom3 and self.atom3 == other.atom1)

    def __lt__(self, other):
        return [self.atom1, self.atom2, self.atom3] < [other.atom1, other.atom2, other.atom3]

    def __gt__(self, other):
        return [self.atom1, self.atom2, self.atom3] > [other.atom1, other.atom2, other.atom3]


class Dihedral():
    def __init__(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

    def __str__(self):
        return '<Dihedral: %s-%s-%s-%s>' \
               % (self.atom1.name, self.atom2.name, self.atom3.name, self.atom4.name)

    def __repr__(self):
        return str(self) + ' at 0x' + str(hex(id(self))[2:].upper())

    def __eq__(self, other):
        return (self.atom1 == other.atom1 and self.atom2 == other.atom2 and
                self.atom3 == other.atom3 and self.atom4 == other.atom4) \
               or (self.atom1 == other.atom4 and self.atom2 == other.atom3 and
                   self.atom3 == other.atom2 and self.atom4 == other.atom1)

    def __lt__(self, other):
        return [self.atom1, self.atom2, self.atom3, self.atom4] \
               < [other.atom1, other.atom2, other.atom3, other.atom4]

    def __gt__(self, other):
        return [self.atom1, self.atom2, self.atom3, self.atom4] \
               > [other.atom1, other.atom2, other.atom3, other.atom4]


class Improper():
    '''
    center atom is the first
    '''

    def __init__(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

    def __str__(self):
        return '<Improper: %s-%s-%s-%s>' \
               % (self.atom1.name, self.atom2.name, self.atom3.name, self.atom4.name)

    def __repr__(self):
        return str(self) + ' at 0x' + str(hex(id(self))[2:].upper())

    def __eq__(self, other):
        if self.atom1 != other.atom1:
            return False

        at12, at13, at14 = sorted([self.atom2, self.atom3, self.atom4])
        at22, at23, at24 = sorted([other.atom2, other.atom3, other.atom4])
        return at12 == at22 and at13 == at23 and at14 == at24

    def __lt__(self, other):
        return [self.atom1, self.atom2, self.atom3, self.atom4] \
               < [other.atom1, other.atom2, other.atom3, other.atom4]

    def __gt__(self, other):
        return [self.atom1, self.atom2, self.atom3, self.atom4] \
               > [other.atom1, other.atom2, other.atom3, other.atom4]
