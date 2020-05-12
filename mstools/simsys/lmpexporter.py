import numpy as np
from .system import System
from ..forcefield.ffterm import *
from ..topology import Topology, Atom, UnitCell, Psf, Bond, Angle, Dihedral, Improper
from ..trajectory import Frame, Trajectory, Gro
from .. import logger


class LammpsExporter():
    def __init__(self):
        pass

    @staticmethod
    def export(system: System, data_out, in_out):
        raise NotImplementedError('This method havn\'t been implemented')
