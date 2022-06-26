from mstk.trajectory.handler import TrjHandler


class CombinedTrj(TrjHandler):
    '''
    Read several trajectory files at the same time.

    It is useful for processing truncated trajectories after restarting simulation.
    All trajectory files should contain the same number of atoms.
    '''

    def __init__(self, files, mode='r'):
        if mode != 'r':
            raise Exception('CombinedTrj is read-only')

        super().__init__()

        self._handlers = []
        for file in files:
            Handler = TrjHandler.get_handler_for_file(file)
            self._handlers.append(Handler(file, mode='r'))

    def get_info(self):
        self._i_frame_handler = []
        self._i_frame_offset = []
        for handler in self._handlers:
            handler.get_info()
            self._i_frame_handler += [handler] * handler.n_frame
            self._i_frame_offset += list(range(handler.n_frame))

        if set([handler.n_atom for handler in self._handlers]).__len__() != 1:
            raise Exception('All trajectories should have same number of atoms')

        self.n_frame = len(self._i_frame_offset)
        self.n_atom = self._handlers[0].n_atom
        return self.n_atom, self.n_frame

    def close(self):
        for handler in self._handlers:
            handler.close()

    def read_frame(self, i_frame, frame):
        handler = self._i_frame_handler[i_frame]
        i = self._i_frame_offset[i_frame]
        handler.read_frame(i, frame)
