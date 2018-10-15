import os
from pathlib import Path
import subprocess
from subprocess import Popen, PIPE

from . import Slurm
from ..errors import JobManagerError


class RemoteSlurm(Slurm):
    def __init__(self, queue, nprocs, ngpu, nproc_request, host, username, remote_dir, port=22, **kwargs):
        super().__init__(queue=queue, nprocs=nprocs, ngpu=ngpu, nprocs_request=nproc_request, **kwargs)

        self.is_remote = True

        self.host = host
        self.port = port
        self.username = username
        self.remote_dir = remote_dir

        # be careful about the port option. -p in ssh, -P in scp
        self._ssh_cmd = f'ssh -o ConnectTimeout=5 -o BatchMode=yes -o StrictHostKeyChecking=no -p {self.port} {self.username}@{self.host}'
        self._scp_cmd = f'scp -o ConnectTimeout=5 -o BatchMode=yes -o StrictHostKeyChecking=no -P {self.port}'
        # require rsync 3.1.0 or above
        self._rsync_cmd = f'rsync --human-readable --recursive --links --perms --executability --times --info=progress2 --rsh="ssh -p {self.port}"'

        if not self.is_working():
            raise JobManagerError('Cannot connect to Slurm @ %s' % host)

    def is_working(self) -> bool:
        cmd = f'{self._ssh_cmd} sinfo --version'
        sp = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
        stdout, stderr = sp.communicate()

        return stdout.decode().startswith('slurm')

    def upload(self, remote_dir=None) -> bool:
        if remote_dir is None:
            remote_dir = self.remote_dir

        print('Upload to', remote_dir)

        mkdir_cmd = f'{self._ssh_cmd} mkdir -p {remote_dir}'
        if subprocess.call(mkdir_cmd, shell=True) != 0:
            return False

        rsync_cmd = f'{self._rsync_cmd} ./ {self.username}@{self.host}:{remote_dir}/'
        if subprocess.call(rsync_cmd, shell=True) != 0:
            return False

        return True

    def download(self, remote_dir=None) -> bool:
        if remote_dir is None:
            remote_dir = self.remote_dir

        print('Download from', remote_dir)
        rsync_cmd = f'{self._rsync_cmd} {self.username}@{self.host}:{remote_dir}/ ./'
        if subprocess.call(rsync_cmd, shell=True) != 0:
            return False

        return True

    def submit(self, sh=None, remote_dir=None, local_dir=None) -> bool:
        if sh is None:
            sh = self.sh

        if remote_dir is None:
            remote_dir = self.remote_dir

        if local_dir is None:
            local_dir = os.getcwd()

        sh_name = Path(sh).name
        remote_sh_name = sh_name + '.remote'
        remote_sh_path = Path(remote_dir, remote_sh_name)

        local_dir_path = Path(local_dir)
        remote_dir_path = Path(remote_dir)

        with open(sh) as f:
            content = f.read()
        with open(remote_sh_name, 'w') as f:
            f.write(content.replace(str(local_dir_path), str(remote_dir_path)))

        scp_cmd = f'{self._scp_cmd} {remote_sh_name} {self.username}@{self.host}:{remote_sh_path}'
        if subprocess.call(scp_cmd, shell=True) != 0:
            return False

        # Be careful. Use "" for several commands over ssh
        sbatch_cmd = f'{self._ssh_cmd} "cd {remote_dir_path}; {self.submit_cmd} {remote_sh_path}"'
        if subprocess.call(sbatch_cmd, shell=True) != 0:
            return False

        return True

    def kill_job(self, name) -> bool:
        job = self.get_job_from_name(name)
        if job is None:
            return False

        cmd = f'{self._ssh_cmd} scancel {job.id}'
        return subprocess.call(cmd, shell=True) == 0

    def get_all_jobs(self):
        # Show all jobs. Then check the user
        cmd = f'{self._ssh_cmd} scontrol show job'
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
