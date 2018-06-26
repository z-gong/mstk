import subprocess
from subprocess import PIPE, Popen
from collections import OrderedDict

from .jobmanager import JobManager
from .node import Node
from .pbsjob import PbsJob
from ..errors import JobManagerError


class Torque(JobManager):
    def __init__(self, queue_list, **kwargs):
        queue = queue_list[0]
        super().__init__(queue=queue[0], nprocs=queue[1], ngpu=queue[2], nprocs_request=queue[3], **kwargs)
        self.sh = '_job_torque.sh'

    def refresh_preferred_queue(self) -> bool:
        return True

        # TODO disable this function
        # if len(self.queue_dict) > 1:
        #     available_queues = self.get_available_queues()
        #     for queue in self.queue_dict.keys():
        #         if queue in available_queues.keys() and available_queues[queue] > 0:
        #             self.queue = queue
        #             self.nprocs = self.queue_dict[queue]
        #             return True
        # 
        # self.queue = list(self.queue_dict.keys())[0]
        # self.nprocs = list(self.queue_dict.values())[0]
        # return False

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

    def submit(self, sh=None):
        if sh is None:
            sh = self.sh
        sp = Popen(['qsub', sh])
        sp.communicate()
        if sp.returncode == 0:
            return True
        else:
            return False

    def kill_job(self, name) -> bool:
        id = self.get_id_from_name(name)
        if id == None:
            return False
        try:
            subprocess.check_call(['qdel', str(id)])
        except:
            raise JobManagerError('Cannot kill job: %s' % name)

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
            job = PbsJob(id=id, name=name, state=state, workdir=workdir, user=user)
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

    def get_nodes(self):
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

    def get_available_queues(self):
        queues = {}
        try:
            nodes = self.get_nodes()
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
