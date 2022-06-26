import os
from io import IOBase


class TrjHandler():
    '''
    Base class of trajectory handlers for various trajectory formats.

    The methods :func:`get_info`, :func:`read_frame` and :func:`write_frame` should be implemented by subclasses.
    The method :func:`close` should also be overriden by subclasses if more works are required more than close the file.
    '''

    _klass_map = {}

    def __init__(self):
        self._file = IOBase()
        self.n_atom = -1
        self.n_frame = -1

    @staticmethod
    def register_format(extension, Handler):
        '''
        Register a handler class for a trajectory format based on the extension name

        Parameters
        ----------
        extension : str
            The extension should start with '.'.
            The upper of lower cases are not distinguished.
        Handler : subclass of TrjHandler
        '''
        TrjHandler._klass_map[extension.lower()] = Handler

    @staticmethod
    def get_handler_for_file(file):
        '''
        Get the appropriate handler class for a trajectory file.

        If a list of string is provided, the CombinedTrj will be returned.
        Otherwise, the handler class will be determined from the extension name of the file.

        Parameters
        ----------
        file : str or list of str
        '''
        from . import CombinedTrj

        if isinstance(file, list):
            Handler = CombinedTrj
        elif isinstance(file, str):
            try:
                Handler = TrjHandler._klass_map[os.path.splitext(file)[-1].lower()]
            except:
                raise Exception('Unknown extension for trajectory file: %s' % file)
        else:
            raise Exception('Only string and list of string are accepted for trajectory file')

        return Handler

    def get_info(self):
        '''
        Get the number of atoms and frames in the trajectory.

        Also record the offset of lines and frames, so that we can read arbitrary frame later.
        It assumes all frames have the same number of atoms.

        Returns
        -------
        n_atom : int
        n_frame : int
        '''
        raise NotImplementedError('Method not implemented')

    def read_frame(self, i_frame, frame):
        '''
        Read a single frame.

        Parameters
        ----------
        i_frame : int
            The index of the frame in the trajectory
        frame : Frame
            The information read from the trajectory will be written into this Frame
        '''
        raise NotImplementedError('Method not implemented')

    def write_frame(self, frame, **kwargs):
        '''
        Write a frame into the trajectory file opened by the handler.

        The handler should be initialized in mode 'w' or 'a'.

        Parameters
        ----------
        frame : Frame
        kwargs : dict
        '''
        raise NotImplementedError('Method not implemented')

    def close(self):
        '''
        Close the handler.
        '''
        self._file.close()
