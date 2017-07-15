import subprocess
from subprocess import Popen
from collections import OrderedDict

from .jobmanager import JobManager


class Slurm(JobManager):
    def __init__(self, queue_dict: OrderedDict):
        # Current only support one queue
        super().__init__(queue=list(queue_dict.keys())[0], nprocs=list(queue_dict.values())[0])
        self.sh = '_job_slurm.sh'

    def generate_sh(self, workdir, commands, name, sh=None):
        if sh is None:
            sh = self.sh
        out = sh[:-2] + 'out'
        err = sh[:-2] + 'err'
        with open(sh, 'w') as f:
            f.write('''#!/bin/bash
#SBATCH --job-name=%(name)s
#SBATCH --partition=%(queue)s
#SBATCH --output=%(out)s
#SBATCH --error=%(err)s
#SBATCH -n %(nprocs)s
#SBATCH --tasks-per-node=%(tasks)s

source /usr/share/Modules/init/bash
unset MODULEPATH
module use /lustre/usr/modulefiles/pi
module purge
module load icc/16.0 impi/5.1 mkl/11.3

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export I_MPI_FABRICS=shm:dapl

cd %(workdir)s
''' % ({'name': name,
        'out': out,
        'err': err,
        'queue': self.queue,
        'nprocs': self.nprocs,
        'tasks': min(self.nprocs, 16),
        'workdir': workdir
        })
                    )
            for cmd in commands:
                f.write(cmd + '\n')

    def submit(self, sh=None):
        if sh is None:
            sh = self.sh
        Popen(['sbatch', sh]).communicate()

    def get_info_from_id(self, id) -> bool:
        # try:
        # output = subprocess.check_output(['qstat', '-f', str(id)])
        # except:
        # return False

        # for line in output.decode().splitlines():
        # if line.strip().startswith('job_state'):
        # state = line.split()[-1]
        # if state in ['R', 'Q']:
        # return True
        # else:
        # return False
        return False

    def get_id_from_name(self, name: str) -> int:
        # try:
        # output = subprocess.check_output(['qstat'])
        # except:
        # raise

        # for line in output.decode().splitlines():
        # if line.find(name) != -1:
        # return int(line.strip().split()[0])

        return None

    def get_info(self, name) -> bool:
        id = self.get_id_from_name(name)
        if id == None:
            return False
        return self.get_info_from_id(id)

    def kill_job(self, name) -> bool:
        # id = self.get_id_from_name(name)
        # if id == None:
        # return False
        # try:
        # subprocess.check_call(['qdel', str(id)])
        # except:
        # raise Exception('Cannot kill job: %s' % name)

        return True
