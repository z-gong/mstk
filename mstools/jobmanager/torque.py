import subprocess
from subprocess import PIPE, Popen

from .jobmanager import JobManager
from .node import Node
from .pbsjob import PbsJob
from ..errors import JobManagerError


class Torque(JobManager):
    '''
    Torque job scheduler.

    Since the development of Torque has stopped since long time ago,
    the support for Torque will be deprecated in the future.

    Parameters
    ----------
    queue : str
        The jobs will be submitted to this queue.
    nprocs : int
        The CPU cores a job can use.
        In most case, it should be equal to :attr:`nprocs_request`.
        If hyper-threading is on, and the simulation software makes good use of hyper-threading, and CGroup is not enabled by the job scheduler,
        it can be two times of :attr:`nprocs_request` for better performance.
    ngpu : int
        The GPU card a job can use.
    nprocs_request : int
        The CPU cores a job will request from job scheduler.
        It must be smaller than the CPU cores on one node.
    env_cmd : str
        The commands for setting up the environment before running real calculations.
        It will be inserted on the top of job scripts.

    Attributes
    ----------
    queue : str
        The jobs will be submitted on this queue.
    nprocs : int
        The CPU cores (or threads if hyper-threading is enabled) a job can use.
    ngpu : int
        The GPU card a job can use.
    nprocs_request : int
        The CPU cores a job will request from job scheduler.
    env_cmd : str
        The commands for setting up the environment before running real calculations.
    sotred_jobs_expire : int
        The lifetime of cached jobs in seconds.
    time : int
        The wall time limit for a job.
    sh : str
        The default name of the job script.
    '''

    #: Whether or not this is a remote job scheduler
    is_remote = False

    def __init__(self, queue, nprocs, ngpu, nprocs_request, **kwargs):
        super().__init__(queue=queue, nprocs=nprocs, ngpu=ngpu, nprocs_request=nprocs_request, **kwargs)
        self.sh = '_job_torque.sh'

    def is_working(self):
        return True

    def generate_sh(self, workdir, commands, name, sh=None, **kwargs):
        if sh is None:
            sh = self.sh
        out = sh[:-2] + 'out'
        err = sh[:-2] + 'err'
        with open(sh, 'w') as f:
            f.write('#!/bin/bash\n'
                    '#PBS -N %(name)s\n'
                    '#PBS -o %(out)s\n'
                    '#PBS -e %(err)s\n'
                    '#PBS -q %(queue)s\n'
                    '#PBS -l walltime=%(time)i:00:00\n'
                    '#PBS -l nodes=1:ppn=%(nprocs_request)s\n\n'
                    '%(env_cmd)s\n\n'
                    'cd %(workdir)s\n\n'
                    % ({'name'          : name,
                        'out'           : out,
                        'err'           : err,
                        'queue'         : self.queue,
                        'time'          : self.time,
                        'nprocs_request': self.nprocs_request,
                        'env_cmd'       : self.env_cmd,
                        'workdir'       : workdir
                        })
                    )
            for cmd in commands:
                f.write(cmd + '\n')

    def submit(self, sh=None, **kwargs):
        if sh is None:
            sh = self.sh
        sp = Popen(['qsub', sh])
        sp.communicate()
        if sp.returncode == 0:
            return True
        else:
            return False

    def kill_job(self, name) -> bool:
        job = self.get_job_from_name(name)
        if job is None:
            return False
        try:
            subprocess.check_call(['qdel', str(job.id)])
        except:
            return False

        return True

    def get_all_jobs(self):
        def get_job_from_str(job_str) -> PbsJob:
            for line in job_str.replace('\n\t', '').splitlines():  # split properties
                line = line.strip()
                if line.startswith('Job Id'):
                    id = int(line.split()[-1])
                    continue
                key, val = line.split(' = ')
                if key == 'Job_Name':
                    name = val
                if key == 'queue':
                    queue = val
                if key == 'Job_Owner':
                    user = val.split('@')[0]  # Job_Owner = username@hostname
                elif key == 'job_state':
                    state_str = val
                    if val == 'Q':
                        state = PbsJob.State.PENDING
                    elif val == 'R':
                        state = PbsJob.State.RUNNING
                    else:
                        state = PbsJob.State.DONE
                elif key == 'init_work_dir':
                    workdir = val
            job = PbsJob(id=id, name=name, state=state, workdir=workdir, user=user, queue=queue)
            job.state_str = state_str
            return job

        # Only show jobs belong to self.username
        cmd = 'qstat -f -u %s' % self.username
        try:
            output = subprocess.check_output(cmd.split())
        except Exception as e:
            raise JobManagerError(str(e))

        jobs = []
        for job_str in output.decode().split('\n\n'):  # split jobs
            if job_str.startswith('Job Id'):
                job = get_job_from_str(job_str)
                # Only show jobs belong to self.username, so this is always True
                if job.user == self.username:
                    jobs.append(job)
        return jobs

    def _get_nodes(self):
        '''
        Deprecated
        '''
        def parse_used_cores(line):
            n_used = 0
            jobs = line.split('=')[-1].strip().split(',')
            for job in jobs:
                cores = job.split('/')[0]
                if '-' in cores:
                    n_used += int(cores.split('-')[1]) - int(cores.split('-')[0]) + 1
                else:
                    n_used += 1
            return n_used

        sp_out = Popen(['pbsnodes'], stdout=PIPE, stderr=PIPE).communicate()
        stdout = sp_out[0].decode()
        stderr = sp_out[1].decode()

        if stderr != '':
            raise Exception('Torque error: pbsnodes failed')

        nodes = []
        for line in stdout.splitlines():
            if not line.startswith(' '):
                name = line.strip()
            elif line.startswith('     state ='):
                state = line.split('=')[1].strip()
            elif line.startswith('     np ='):
                np = int(line.split('=')[1].strip())
                n_used_cores = 0
            elif line.startswith('     properties = '):
                queue = line.split('=')[1].strip()
            elif line.startswith('     jobs = '):
                n_used_cores = parse_used_cores(line)
            elif line.startswith('     status = '):
                n_free_cores = np - n_used_cores
                nodes.append(Node(name, queue, np, n_free_cores, state))

        return nodes

    def _get_available_queues(self):
        '''
        Deprecated
        '''
        queues = {}
        try:
            nodes = self._get_nodes()
        except Exception as e:
            print(str(e))
        else:
            for node in nodes:
                if node.queue not in queues.keys():
                    queues[node.queue] = 0
                if node.state == 'free' and node.queue in self.queue_dict.keys() \
                        and node.n_free_cores >= self.queue_dict[node.queue]:
                    queues[node.queue] += 1

        return queues
