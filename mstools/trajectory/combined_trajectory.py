from .trajectory import Trajectory


class CombinedTrajectory(Trajectory):
    def __init__(self, files, mode='r'):
        if mode != 'r':
            raise Exception('CombinedTrajectory is read-only')

        super().__init__()

        self._trajectories = []
        for file in files:
            self._trajectories.append(Trajectory.open(file, mode='r'))

        if set([trj.n_atom for trj in self._trajectories]).__len__() != 1:
            raise Exception('All trajectories should have same number of atoms')

        self._i_frame_in_trj = []
        self._i_frame_offset = []
        for trj in self._trajectories:
            self._i_frame_in_trj += [trj] * trj.n_frame
            self._i_frame_offset += list(range(trj.n_frame))

        self.n_frame = len(self._i_frame_offset)
        self.n_atom = self._trajectories[0].n_atom

        self._mode = mode
        self._opened = True

    def close(self):
        for trj in self._trajectories:
            trj.close()
        self._opened = False

    def _read_frame(self, i_frame: int, frame):
        trj = self._i_frame_in_trj[i_frame]
        i = self._i_frame_offset[i_frame]
        return trj._read_frame(i, frame)
