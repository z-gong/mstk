import numpy as np
from mstk import logger
from mstk.topology import UnitCell


class Frame:
    '''
    A frame records the positions and unit cell (and optionally velocities, charges, current step and current time)
    of the system at one specific time step.

    Velocity is stored for analyzing dynamic properties with auto correlation function.
    Charge is stored for fluctuating charge method.

    Parameters
    ----------
    n_atom : int
        The number of atoms in the frame.

    Attributes
    ----------
    step : int
        The current step. -1 means unknown.
    time : float
        The current simulation time. -1 means unknown.
    cell : UnitCell
        The unit cell.
    positions : array_like
        The positions of all atoms in the frame.
        It's a np.ndarry of shape (n_atom, 3)
    has_velocity : bool
        Whether or not there's velocity information for all atoms
    velocities : array_like
        The velocities of all atoms in the frame.
        It's a np.ndarry of shape (n_atom, 3)
    has_charge : bool
        Whether or not there's charge information for all atoms
    charges : array_like
        The charge of all atoms in the frame.
        It's a np.ndarry of shape (n_atom,)
    '''

    def __init__(self, n_atom):
        self.step = -1  # -1 means step information not found
        self.time = -1  # in ps. -1 means time information not found
        self.cell = UnitCell()
        self.n_atom = n_atom
        self.has_velocity = False
        self.has_charge = False
        self._positions = np.zeros((n_atom, 3), dtype=float)
        self._velocities = np.zeros((n_atom, 3), dtype=float)
        self._charges = np.zeros(n_atom, dtype=float)  # for fluctuating charge simulations

    def reset(self):
        '''
        Remove step, time and cell information, and reset has_velocity, has_charge to False.
        '''
        self.step = -1
        self.time = -1
        self.cell.set_box(None)
        self.has_velocity = False
        self.has_charge = False

    def resize(self, n_atom):
        '''
        Resize the number of atoms recorded in the frame.

        Usually it is not necessary to call this method explicitly.
        If n_atom is smaller than previous value, the information of the end atoms will be trimmed.
        If n_atom is larger than previous value, the information for the new atoms will be random.

        Parameters
        ----------
        n_atom : int
            The new number of atoms
        '''
        self.n_atom = n_atom
        self._positions.resize((n_atom, 3), refcheck=False)
        self._velocities.resize((n_atom, 3), refcheck=False)
        self._charges.resize((n_atom, 3), refcheck=False)

    @property
    def positions(self):
        '''
        The positions of all atoms in this frame

        :setter: Set the positions of all atoms

        Returns
        -------
        positions : array_like
            The positions is a numpy array of shape (n_atom, 3)
        '''
        return self._positions

    @positions.setter
    def positions(self, value):
        if not isinstance(value, (list, tuple, np.ndarray)) or len(value) != self.n_atom \
                or not isinstance(value[0], (list, tuple, np.ndarray)) or len(value[0]) != 3:
            raise ValueError(f'positions should have the shape of ({self.n_atom}, 3)')
        self._positions[:] = value

    @property
    def velocities(self):
        '''
        The velocities of all atoms in this frame

        :setter: Set the velocities of all atoms

        Returns
        -------
        velocities : array_like
            The velocities is a numpy array of shape (n_atom, 3)
        '''
        return self._velocities

    @velocities.setter
    def velocities(self, value):
        if not isinstance(value, (list, tuple, np.ndarray)) or len(value) != self.n_atom \
                or not isinstance(value[0], (list, tuple, np.ndarray)) or len(value[0]) != 3:
            raise ValueError(f'velocities should have the shape of ({self.n_atom}, 3)')
        self._velocities[:] = value

    @property
    def charges(self):
        '''
        The charges of all atoms in this frame

        :setter: Set the charges of all atoms

        Returns
        -------
        charges : array_like
            The charges is a numpy array of shape (n_atom,)
        '''
        return self._charges

    @charges.setter
    def charges(self, value):
        if not isinstance(value, (list, tuple, np.ndarray)) or len(value) != self.n_atom:
            raise ValueError(f'charges should have {self.n_atom} elements')
        self._charges[:] = value
