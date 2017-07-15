import subprocess
from subprocess import PIPE, Popen
from collections import OrderedDict

from .jobmanager import JobManager
from .node import Node


class Torque(JobManager):
    def __init__(self, queue_dict: OrderedDict):
        super().__init__(list(queue_dict.keys())[0], list(queue_dict.values())[0])
        self.queue_dict = queue_dict
        self.sh = '_job_torque.sh'

    def refresh_preferred_queue(self) -> bool:
        if len(self.queue_dict) > 1:
            available_queues = self.get_available_queues()
            for queue in self.queue_dict.keys():
                if queue in available_queues.keys() and available_queues[queue] > 0:
                    self.queue = queue
                    self.nprocs = self.queue_dict[queue]
                    return True

        self.queue = list(self.queue_dict.keys())[0]
        self.nprocs = list(self.queue_dict.values())[0]
        return False

    def generate_sh(self, workdir, commands, name, sh=None):
        if sh is None:
            sh = self.sh
        out = sh[:-2] + 'out'
        err = sh[:-2] + 'err'
        with open(sh, 'w') as f:
            f.write('#!/bin/sh\n'
                    '#PBS -N %(name)s\n'
                    '#PBS -o %(out)s\n'
                    '#PBS -e %(err)s\n'
                    '#PBS -q %(queue)s\n'
                    '#PBS -l nodes=1:ppn=%(nprocs)s\n\n'
                    'cd %(workdir)s\n'
                    % ({'name': name,
                        'out': out,
                        'err': err,
                        'queue': self.queue,
                        'nprocs': self.nprocs,
                        'workdir': workdir
                        })
                    )
            for cmd in commands:
                f.write(cmd + '\n')

    def submit(self, sh=None):
        if sh is None:
            sh = self.sh
        Popen(['qsub', sh]).communicate()

    def get_info_from_id(self, id) -> bool:
        try:
            output = subprocess.check_output(['qstat', '-f', str(id)])
        except:
            return False

        for line in output.decode().splitlines():
            if line.strip().startswith('job_state'):
                state = line.split()[-1]
                if state in ['R', 'Q']:
                    return True
                else:
                    return False
        return False

    def get_id_from_name(self, name: str) -> int:
        try:
            output = subprocess.check_output(['qstat'])
        except:
            raise

        for line in output.decode().splitlines():
            if line.find(name) != -1:
                return int(line.strip().split()[0])

        return None

    def get_info(self, name) -> bool:
        id = self.get_id_from_name(name)
        if id == None:
            return False
        return self.get_info_from_id(id)

    def kill_job(self, name) -> bool:
        id = self.get_id_from_name(name)
        if id == None:
            return False
        try:
            subprocess.check_call(['qdel', str(id)])
        except:
            raise Exception('Cannot kill job: %s' % name)

        return True

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
