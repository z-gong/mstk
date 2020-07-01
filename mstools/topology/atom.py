import numpy as np


class Atom():
    '''
    An atom is a particle in simulation.

    It can be an real atom, a Drude particle, a virtual site or a coarse-grained bead.
    '''

    def __init__(self, name='UNK'):
        '''
        Initialize a atom with name.

        The full information will be assigned later manually or by force field.

        Parameters
        ----------
        name : str
        '''
        #: name of this atom, not necessarily unique
        self.name = name
        #: id in topology. -1 means information haven\'t been updated by topology
        self.id = -1
        #: atom type in force field. Required for assigning force field parameters
        self.type = ''
        #: atomic symbol. Mainly used for output topology or trajectory
        self.symbol = 'UNK'
        #: atomic mass. Required for output simulation files
        self.mass = 0.
        #: atomic mass. Required for output simulation files
        self.charge = 0.
        #: isotropic polarizability. Required for output Drude polarizable simulation files
        self.alpha = 0.
        #: thole screening for induced dipole. Required for output Drude polarizable simulation files
        self.thole = 0.
        #: formal charge calculate from valence bond theory. Optionally required by typing engine
        self.formal_charge = 0.
        #: whether or not this is a Drude particle for polarizable model
        self.is_drude = False
        #: the instance of subclass of VirtualSite if this atom is a virtual site
        self.virtual_site = None
        #: whether or not this atom contains position
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
        '''
        The molecule this atom belongs to

        Returns
        -------
        molecule : Molecule
        '''
        return self._molecule

    @property
    def id_in_molecule(self):
        '''
        The index of this atom in molecule.

        Returns
        -------
        index : int
            -1 means information haven\'t been assigned by molecule or topology
        '''
        return self._id_in_molecule

    @property
    def bonds(self):
        '''
        All the bonds which contains this atom.

        Bonds of Drude dipoles are included if exist.

        Returns
        -------
        bonds : list of Bond
        '''
        return self._bonds

    @property
    def bond_partners(self):
        '''
        All the bond partners.

        Drude particles are included if exist.

        Returns
        -------
        partners : list of Atom
        '''
        return [bond.atom2 if bond.atom1 == self else bond.atom1 for bond in self._bonds]

    @property
    def position(self):
        '''
        The position of this atom

        :setter: Set the position of this atom

        Returns
        -------
        position : array_like
        '''
        return self._position

    @position.setter
    def position(self, value):
        if not isinstance(value, (list, tuple, np.ndarray)) or len(value) != 3:
            raise ValueError('position should has three elements')
        self._position[:] = value
        self.has_position = True
