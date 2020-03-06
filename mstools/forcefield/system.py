from mstools.topology.topology import Topology
from .parameterset import ParameterSet
from ..trajectory import Frame


class System():
    def __init__(self, topology: Topology, params: ParameterSet, frame: Frame):
        self._topology = topology
        self._params = params
        self._frame = frame

    def write_lmp(self):
        pass

    def write_gmx(self):
        pass

    def write_omm(self):
        pass
