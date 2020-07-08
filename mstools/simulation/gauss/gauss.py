import os
import shutil

from ...errors import GmxError
from ...wrapper import Gauss
from ...utils import get_last_line
from ...util import DocstringMeta


class GaussSimulation(metaclass=DocstringMeta):
    '''
    Base class of predefined simulation protocols with Gaussian.

    The following methods should be implemented/overridden by its subclasses:

    * :func:`set_system`
    * :func:`prepare`
    * :func:`run`
    * :func:`check_finished`
    * :func:`analyze`
    * :func:`clean`

    GaussSimulation should not be constructed directly. Use its subclasses instead.
    '''

    def __init__(self, gauss=None, gauss_bin=None, gauss_scrdir=None, jobmanager=None, **kwargs):
        if gauss is not None:
            self.gauss  = gauss
        elif gauss_bin is not None:
            self.gauss = Gauss(gauss_bin=gauss_bin, scrdir=gauss_scrdir)
        else:
            raise Exception('argument gauss missing')

        if jobmanager is None:
            raise Exception('argument jobmanager missing')
        self.jobmanager = jobmanager

        self.logs = []  # used for checking whether the job is successfully finished

    def set_system(self, smiles_list, n_mol_list):
        '''
        Specify the components of the simulation system.

        Parameters
        ----------
        smiles_list : list of str
            List of SMILES of the molecules to be simulated
        n_mol_list : list of int, optional
            List of numbers of molecules
        '''
        if type(smiles_list) != list:
            raise Exception('smiles_list should be list')
        self.smiles_list = smiles_list[:]
        self.n_mol_list = n_mol_list[:]

    def prepare(self, **kwargs):
        '''
        Generate input files for simulation engines, and return the commands for running the simulation.

        This method should be implemented by subclasses.

        Returns
        -------
        commands : list of str
            The commands for running the simulation.
        '''
        raise NotImplementedError('Method not implemented')

    def run(self):
        '''
        Submit the job to PBS scheduler
        '''
        self.jobmanager.submit()

    def check_finished(self):
        '''
        Check if the job is successfully finished.

        This method should be implemented by subclasses.
        '''
        raise NotImplementedError('Method not implemented')

    def analyze(self):
        '''
        Analyze simulation data and obtain raw results like energy, density...

        This method should be implemented by subclasses.
        '''
        raise NotImplementedError('Method not implemented')

    def clean(self):
        '''
        Delete intermediate simulation data.

        This method should be implemented by subclasses.
        '''
        raise NotImplementedError('Method not implemented')
