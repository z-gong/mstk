import math
import subprocess
from subprocess import Popen, PIPE

from .pbsjob import PbsJob
from .jobmanager import JobManager
from ..errors import JobManagerError


class Slurm(JobManager):
    '''
    Slurm job scheduler.

    '''
    def __init__(self, queue, nprocs, ngpu, nprocs_request, **kwargs):
        super().__init__(queue=queue, nprocs=nprocs, ngpu=ngpu, nprocs_request=nprocs_request, **kwargs)
        self.sh = '_job_slurm.sh'
        self.submit_cmd = 'sbatch'

    def is_working(self) -> bool:
        cmd = 'sinfo --version'
        sp = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
        stdout, stderr = sp.communicate()
        return stdout.decode().startswith('slurm')

    def replace_mpirun_srun(self, commands) -> (int, [str]):
        n_mpi = 1
        cmds_replaced = []
        for cmd in commands:
            if cmd.startswith('mpirun'):
                n_mpi = int(cmd.split()[2])
                cmd_srun = 'srun -n %i ' % n_mpi + ' '.join(cmd.split()[3:])
                cmds_replaced.append(cmd_srun)
            else:
                cmds_replaced.append(cmd)
        return n_mpi, cmds_replaced

    def generate_sh(self, workdir, commands, name, sh=None, n_tasks=None, **kwargs):
        '''
        If do not set n_tasks, default a whole node is used (which means n_tasks equal to self.nprocs_request)

        :param workdir:
        :param commands:
        :param name:
        :param sh:
        :param n_tasks:
        :param kwargs:
        :return:
        '''
        if sh is None:
            sh = self.sh
        out = sh[:-2] + 'out'
        err = sh[:-2] + 'err'

        n_mpi, srun_commands = self.replace_mpirun_srun(commands)

        if n_tasks is None:
            n_node = 1
            n_tasks = self.nprocs_request
        else:
            n_node = int(math.ceil(n_tasks / self.nprocs_request))

        if self.ngpu > 0:
            gpu_cmd = '#SBATCH --gres=gpu:%i\n' % self.ngpu * n_node
        else:
            gpu_cmd = ''

        with open(sh, 'w') as f:
            f.write('#!/bin/bash\n'
                    '#SBATCH -D %(workdir)s\n'
                    '#SBATCH -J %(name)s\n'
                    '#SBATCH -o %(out)s\n'
                    '#SBATCH -e %(err)s\n'
                    '#SBATCH -p %(queue)s\n'
                    '#SBATCH --time=%(time)i:00:00\n'
                    '#SBATCH --nodes=%(n_node)i\n'
                    '#SBATCH --ntasks=%(n_tasks)i\n'
                    '%(gpu_cmd)s'
                    '\n'
                    '%(env_cmd)s\n\n'
                    % ({'name'   : name,
                        'out'    : out,
                        'err'    : err,
                        'queue'  : self.queue,
                        'time'   : self.time,
                        'n_node' : n_node,
                        'n_tasks': n_tasks,
                        'gpu_cmd': gpu_cmd,
                        'env_cmd': self.env_cmd,
                        'workdir': workdir
                        })
                    )
            for cmd in srun_commands:
                f.write(cmd + '\n')

    def submit(self, sh=None, **kwargs) -> bool:
        if sh is None:
            sh = self.sh
        cmd = self.submit_cmd + ' ' + sh
        return subprocess.call(cmd.split()) == 0

    def kill_job(self, name) -> bool:
        job = self.get_job_from_name(name)
        if job is None:
            return False

        cmd = f'scancel {job.id}'
        return subprocess.call(cmd.split()) == 0

    def get_all_jobs(self):
        # Show all jobs. Then check the user
        cmd = 'scontrol show job'
        sp = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
        stdout, stderr = sp.communicate()
        if sp.returncode != 0:
            print(stderr.decode())
            return []

        jobs = []
        for job_str in stdout.decode().split('\n\n'):  # split jobs
            if job_str.startswith('JobId'):
                job = self.get_job_from_str(job_str)
                # Show all jobs. Then check the user
                if job.user == self.username:
                    jobs.append(job)
        return jobs

    def get_job_from_str(self, job_str) -> PbsJob:
        workdir = None
        for line in job_str.split():  # split properties
            try:
                key, val = line.split('=')[0:2]
            except:
                continue
            if key == 'JobId':
                id = int(val)
            elif key == 'UserId':
                user = val.split('(')[0]  # UserId=username(uid)
            elif key == 'JobName' or key == 'Name':
                name = val
            elif key == 'Partition':
                queue = val
            elif key == 'JobState':
                state_str = val
                if val in ('PENDING', 'RESV_DEL_HOLD'):
                    state = PbsJob.State.PENDING
                elif val in ('CONFIGURING', 'RUNNING', 'COMPLETING', 'STOPPED', 'SUSPENDED'):
                    state = PbsJob.State.RUNNING
                else:
                    state = PbsJob.State.DONE
            elif key == 'WorkDir':
                workdir = val
        job = PbsJob(id=id, name=name, state=state, workdir=workdir, user=user, queue=queue)
        job.state_str = state_str
        return job
