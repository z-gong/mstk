from ..topology import Topology
from ..trajectory import Frame, Trajectory
from .parameterset import ParameterSet


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

    def write_gro(self, file):
        from ..trajectory import Gro
        gro = Gro(file, 'w')
        gro.write_frame(self._topology, self._frame)
        gro.close()

    def write_psf(self, file):
        from ..topology import Psf
        psf = Psf(file, 'w')
        psf.is_drude = self._topology.is_drude
        psf.init_from_molecules(self._topology.molecules)
        psf.write()
