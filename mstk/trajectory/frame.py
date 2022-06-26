import numpy as np
from mstk import logger
from mstk.topology import UnitCell


class Frame():
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
        self.positions = np.zeros((n_atom, 3), dtype=np.float32)
        self.has_velocity = False
        self.velocities = np.zeros((n_atom, 3), dtype=np.float32)
        self.has_charge = False
        self.charges = np.zeros(n_atom, dtype=np.float32)  # for fluctuating charge simulations

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
        If n_atom is larger than previous value, the information for the new atoms will be random .

        Parameters
        ----------
        n_atom : int
            The new number of atoms
        '''
        if n_atom < len(self.positions):
            logger.warning('n_atom is smaller than original. '
                           'The positions, velocities and charges at the end will be lost')
        self.positions.resize((n_atom, 3), refcheck=False)
        self.velocities.resize((n_atom, 3), refcheck=False)
        self.charges.resize((n_atom, 3), refcheck=False)
