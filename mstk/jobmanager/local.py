from .jobmanager import JobManager
from ..errors import JobManagerError


class Local(JobManager):
    '''
    A fake job scheduler that run jobs on the local machine.

    It is designed for generate simulation files locally for debugging purpose.

    Parameters
    ----------
    nprocs : int
    ngpu : int
    env_cmd : str

    Attributes
    ----------
    nprocs : int
    ngpu : int
    nprocs_request : int
    env_cmd : str
    username : str
    sh : str
    '''

    #: Whether or not this is a remote job scheduler
    is_remote = False

    def __init__(self, nprocs, ngpu, **kwargs):
        super().__init__(nprocs=nprocs, ngpu=ngpu, **kwargs)
        self.sh = '_job_local.sh'

    def is_working(self):
        return True

    def generate_sh(self, workdir, commands, name=None, sh=None, **kwargs):
        if sh is None:
            sh = self.sh
        with open(sh, 'w') as f:
            f.write('#!/bin/sh\n\n'
                    'cd %(workdir)s\n'
                    % ({'workdir': workdir})
                    )
            for cmd in commands:
                f.write(cmd + '\n')

    def submit(self, sh=None):
        raise JobManagerError('Not supported on localhost')

    def is_running(self, name):
        raise JobManagerError('Not supported on localhost')

    def kill_job(self, name):
        raise JobManagerError('Not supported on localhost')

    def get_all_jobs(self):
        return []
