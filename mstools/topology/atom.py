import numpy as np


class Atom():
    def __init__(self, name='UNK'):
        self.id = -1  # id in topology
        self.name = name
        self.type = ''
        self.symbol = 'UNK'
        self.mass = 0.
        self.charge = 0.
        self.alpha = 0.  # isotropic polarizability
        self.thole = 0.  # thole screening for induced dipole
        self.formal_charge = 0.
        self.is_drude = False  # this is a Drude particle for polarizable model
        self.virtual_site = None  # this is a virtual site
        self.has_position = False
        self._position = np.zeros(3, dtype=np.float32)
        self._molecule = None
        self._bonds = []
        self._id_in_molecule = -1

    def __repr__(self):
        return f'<Atom: {self.name} {self.id} {self.type}>'

    def __lt__(self, other):
        return self.id < other.id

    def __gt__(self, other):
        return self.id > other.id

    def __deepcopy__(self, memodict={}):
        '''
        VirtualSite information is not copied
        Because it makes no sense to have a single VirtualSite atom without the knowledge of parent atoms
        '''
        atom = Atom()
        atom.id = self.id
        atom.name = self.name
        atom.type = self.type
        atom.symbol = self.symbol
        atom.mass = self.mass
        atom.charge = self.charge
        atom.formal_charge = self.formal_charge
        atom.is_drude = self.is_drude
        atom.has_position = self.has_position
        atom._position[:] = self._position[:]
        atom.alpha = self.alpha
        atom.thole = self.thole
        return atom

    @property
    def molecule(self):
        return self._molecule

    @property
    def id_in_molecule(self):
        return self._id_in_molecule

    @property
    def bonds(self):
        return self._bonds

    @property
    def bond_partners(self):
        return [bond.atom2 if bond.atom1 == self else bond.atom1 for bond in self._bonds]

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, value):
        if not isinstance(value, (list, tuple, np.ndarray)) or len(value) != 3:
            raise ValueError('position should has three elements')
        self._position[:] = value
        self.has_position = True
