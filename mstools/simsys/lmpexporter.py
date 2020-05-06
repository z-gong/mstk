import numpy as np
import warnings
from .system import System
from ..forcefield.ffterm import *
from ..topology import Topology, Atom, UnitCell, Psf, Bond, Angle, Dihedral, Improper
from ..trajectory import Frame, Trajectory, Gro


class LammpsExporter():
    def __init__(self):
        pass

    @staticmethod
    def export(system: System, data_out, in_out):
        raise NotImplementedError('This method havn\'t been implemented')
