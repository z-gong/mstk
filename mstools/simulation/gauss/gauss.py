import os
import shutil

from ..simulation import Simulation
from ...errors import GmxError
from ...wrapper import Gauss
from ...utils import get_last_line


class GaussSimulation(Simulation):
    '''
    Base class of predefined simulation protocols with Gaussian.

    Parameters
    ----------
    gauss_bin : str
        The binary of Gaussian
    gauss_scrdir : str
        The scratch directory to be used by Gaussian calculation
    jobmanager : subclass of JobManager
    '''
    def __init__(self, gauss_bin, gauss_scrdir=None, jobmanager=None, **kwargs):
        super().__init__(jobmanager=jobmanager, **kwargs)
        self.gauss = Gauss(gauss_bin=gauss_bin, scrdir=gauss_scrdir)
        self.logs = []  # used for checking whether the job is successfully finished

