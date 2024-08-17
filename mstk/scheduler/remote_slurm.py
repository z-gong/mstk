import os
import subprocess
from subprocess import Popen, PIPE
from pathlib import Path
from mstk import logger
from mstk.errors import SchedulerError
from . import Slurm


class RemoteSlurm(Slurm):
    '''
    Slurm job scheduler running on a remote machine.

    Parameters
    ----------
    host : str
        The IP address of the remote host that is running Slurm
    port : int
        The SSH port for logging in the remote host
    username : str
        The username for logging in the remote host
    remote_dir : str
        The default directory to use on the remote host for running calculation

    Attributes
    ----------
    host : str
        The IP address of the remote host that is running Slurm
    port : int
        The SSH port for logging in the remote host
    username : str
        The username for logging in the remote host
    remote_dir : str
        The default directory to use on the remote host for running calculation
    sh : str
        The default name of the job script
    job_parameter : JobParameter
        The default parameters for submitting a job
    submit_cmd : str
        The command for submitting the job script.
        If is `sbatch` by default. But extra argument can be provided, e.g. `sbatch --qos=debug`.
    cached_jobs_expire : int
        The lifetime of cached jobs in seconds.
    '''

    #: Whether this is a remote job scheduler
    is_remote = True

    def __init__(self, host, username, remote_dir, port=22):
        super().__init__()

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
            raise SchedulerError('Cannot connect to Slurm @ %s' % host)

    def is_working(self) -> bool:
        '''
        Check whether Slurm is working normally on the remote machine.

        It calls `sinfo --version` and check the output.

        Returns
        -------
        is : bool
        '''
        logger.info(f'Connecting to Slurm on {self.host}')
        cmd = f'{self._ssh_cmd} sinfo --version'
        sp = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
        stdout, stderr = sp.communicate()

        return stdout.decode().startswith('slurm')

    def upload(self, local_dir=None, remote_dir=None):
        '''
        Upload all the files in current local directory to remote directory.

        Parameters
        ----------
        local_dir : dir, optional
            If not set, will use the current dir.
        remote_dir : dir, optional
            If not set, will use the default :attr:`remote_dir`.

        Returns
        -------
        successful : bool
            Whether the upload is successful
        '''

        local_dir = local_dir or os.getcwd()
        remote_dir = remote_dir or self.remote_dir

        logger.info(f'Uploading to {remote_dir}')

        mkdir_cmd = f'{self._ssh_cmd} mkdir -p {remote_dir}'
        if subprocess.call(mkdir_cmd, shell=True) != 0:
            return False

        rsync_cmd = f'{self._rsync_cmd} {local_dir}/ {self.username}@{self.host}:{remote_dir}/'
        if subprocess.call(rsync_cmd, shell=True) != 0:
            return False

        return True

    def download(self, remote_dir=None, local_dir=None) -> bool:
        '''
        Upload all the files in remote directory to current local directory.

        Parameters
        ----------
        remote_dir : dir, optional
            If not set, will use the default :attr:`remote_dir`.
        local_dir : dir, optional
            If not set, will use the current dir.

        Returns
        -------
        successful : bool
            Whether the download is successful
        '''
        remote_dir = remote_dir or self.remote_dir
        local_dir = local_dir or os.getcwd()

        logger.info(f'Downloading from {remote_dir}')
        rsync_cmd = f'{self._rsync_cmd} {self.username}@{self.host}:{remote_dir}/ {local_dir}/'
        if subprocess.call(rsync_cmd, shell=True) != 0:
            return False

        return True

    def submit(self, sh=None, remote_dir=None):
        '''
        Submit a job script to the Slurm scheduler on the remote machine.

        Parameters
        ----------
        sh : str
            The job script to be submitted.
        remote_dir : str
            The directory to submit the script on the remote machine.

        Returns
        -------
        id : int
        '''
        sh = sh or self.sh
        sh_name = Path(sh).name

        remote_dir_path = Path(remote_dir or self.remote_dir).absolute()
        remote_sh_path = remote_dir_path / sh_name

        logger.info(f'Uploading slurm script to {remote_dir_path}')

        mkdir_cmd = f'{self._ssh_cmd} "mkdir -p {remote_dir_path}"'
        if subprocess.call(mkdir_cmd, shell=True) != 0:
            return False

        scp_cmd = f'{self._scp_cmd} {sh} {self.username}@{self.host}:{remote_sh_path}'
        if subprocess.call(scp_cmd, shell=True) != 0:
            return False

        # Be careful. Use "" for several commands over ssh
        sbatch_cmd = f'{self._ssh_cmd} "cd {remote_dir_path}; {self.submit_cmd} {remote_sh_path}"'
        p = Popen(sbatch_cmd, shell=True, stdout=PIPE, stderr=PIPE)
        b_out, b_err = p.communicate()
        out, err = b_out.strip().decode(), b_err.strip().decode()
        if p.returncode == 0 and out.startswith('Submitted batch job'):
            logger.info(f'{out}')
            return int(out.split()[-1])
        else:
            if out:
                logger.error(f'{out}')
            if err:
                logger.error(f'{err}')
            if not out and not err:
                logger.error(f'sbatch failed')
            return -1

    def kill_job(self, name) -> bool:
        job = self.get_job_from_name(name)
        if job is None:
            return False

        cmd = f'{self._ssh_cmd} scancel {job.id}'
        return subprocess.call(cmd, shell=True) == 0

    def _get_jobs(self):
        # Show all jobs. Then check the user
        logger.info(f'Querying jobs on remote Slurm...')

        cmd = f'{self._ssh_cmd} scontrol show job'
        sp = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
        stdout, stderr = sp.communicate()
        if sp.returncode != 0:
            logger.error(stderr.decode())
            return []

        jobs = []
        for job_str in stdout.decode().split('\n\n'):  # split jobs
            if job_str.startswith('JobId'):
                job = self._get_job_from_str(job_str)
                # Show all jobs. Then check the user
                if job.user == self.username:
                    jobs.append(job)
        return jobs
