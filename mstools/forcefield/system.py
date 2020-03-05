from .topology import Topology
from .parameter import ParameterSet
from ..trajectory import Frame


class System():
    def __init__(self, topology: Topology, parameter: ParameterSet, frame: Frame):
        self._topology = topology
        self._parameter = parameter
        self._frame = frame

    def write_lmp(self):
        pass

    def write_gmx(self):
        pass

    def write_omm(self):
        pass
