import numpy as np


class Atom():
    '''
    An atom is a particle in simulation.

    It can be an real atom, a Drude particle, a virtual site or a coarse-grained bead.

    Parameters
    ----------
    name : str

    Attributes
    ----------
    id : int
        Index of this atom in topology. -1 means information haven\'t been updated by topology
    id_in_mol : int
        Index of this atom in molecule. -1 means information haven\'t been updated by topology
    name : str
        Name of this atom, not necessarily unique
    type : str
        Atom type in force field. Required for assigning force field parameters
    symbol : str
        Atomic symbol. Mainly used for output topology or trajectory
    mass : float
        Mass assigned. Required for output simulation files
    charge : float
        Charge assigned. Required for output simulation files
    alpha : float
        Isotropic polarizability. Required for output Drude polarizable simulation files
        This property is set for parent atom, not Drude particle
    thole : float
        Thole screening for induced dipole. Required for output Drude polarizable simulation files
        This property is set for parent atom, not Drude particle
    formal_charge : int
        Formal charge calculated from valence bond theory. Optionally required by typing engine
    is_drude : bool
        Whether or not this is a Drude particle for polarizable model
    virtual_site : None or subclass of VirtualSite
        None if this atom is not a virtual site.
        Otherwise the instance of subclass of VirtualSite
    has_position : bool
        Whether or not the position of this atom is set
    '''

    def __init__(self, name='UNK'):
        self.id = -1
        self.id_in_mol = -1
        self.name = name
        self.type = 'UNK'
        self.symbol = 'UNK'
        self.mass = 0.
        self.charge = 0.
        self.alpha = 0.  # this property is set for parent atom, not Drude particle
        self.thole = 0.  # this property is set for parent atom, not Drude particle
        self.formal_charge = 0
        self.is_drude = False
        self.virtual_site = None
        self.has_position = False

        self._position = np.zeros(3, dtype=float)
        self._molecule = None
        self._bonds = []
        self._residue = None

    def __repr__(self):
        return f'<Atom: {self.name} {self.id} {self.type}>'

    def __lt__(self, other):
        return self.id < other.id

    def __gt__(self, other):
        return self.id > other.id

    def __deepcopy__(self, memodict={}):
        '''
        id, id_in_mol, virtual_site, residue are not copied, because these information depends on other atoms.
        '''
        atom = Atom()
        atom.name = self.name
        atom.type = self.type
        atom.symbol = self.symbol
        atom.mass = self.mass
        atom.charge = self.charge
        atom.alpha = self.alpha
        atom.thole = self.thole
        atom.formal_charge = self.formal_charge
        atom.is_drude = self.is_drude
        atom.has_position = self.has_position
        atom._position[:] = self._position[:]
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
    def residue(self):
        '''
        The residue this atom belongs to

        Returns
        -------
        residue : Residue
        '''
        return self._residue

    @property
    def bonds(self):
        '''
        All the bonds involving this atom.

        Bonds of Drude dipoles are included if exist.

        Returns
        -------
        bonds : list of Bond
        '''
        return self._bonds

    @property
    def bond_partners(self):
        '''
        All the bond partners of this atom.

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
            The position is a numpy array of shape (3,)
        '''
        return self._position

    @position.setter
    def position(self, value):
        if not isinstance(value, (list, tuple, np.ndarray)) or len(value) != 3:
            raise ValueError('position should has three elements')
        self._position[:] = value
        self.has_position = True

    @property
    def rdatom(self):
        '''
        The `rdkit.Chem.Atom` object associated with this atom.

        It exists only if this atom belongs to a molecule.

        Returns
        -------
        rdbond : rdkit.Chem.Atom
        '''
        if self.molecule is None:
            return None

        return self.molecule.rdmol.GetAtomWithIdx(self.id_in_mol)
