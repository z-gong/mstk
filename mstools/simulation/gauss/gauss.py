import os
import shutil

from ..simulation import Simulation
from ...errors import GmxError
from ...wrapper import Gauss
from ...utils import get_last_line


class GaussSimulation(Simulation):
    '''
    Base class of predefined simulation protocols with Gaussian.

    GaussSimulation should not be constructed directly. Use its subclasses instead.
    '''
    def __init__(self, gauss_bin, gauss_scrdir=None, jobmanager=None, **kwargs):
        super().__init__(jobmanager=jobmanager, **kwargs)
        self.gauss = Gauss(gauss_bin=gauss_bin, scrdir=gauss_scrdir)
        self.logs = []  # used for checking whether the job is successfully finished

