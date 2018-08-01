from .jobmanager import JobManager
from ..errors import JobManagerError


class Local(JobManager):
    def __init__(self, nprocs, ngpu, **kwargs):
        super().__init__(nprocs=nprocs, ngpu=ngpu, **kwargs)
        self.sh = '_job_local.sh'

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
        print('Localhost is only for test')
        return True

    def is_running(self, name):
        raise JobManagerError('Not supported on localhost')

    def kill_job(self, name):
        raise JobManagerError('Not supported on localhost')

    def get_all_jobs(self):
        return []
