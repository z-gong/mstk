from .job import Job


class JobManager:
    def __init__(self, queue=None, nprocs=1, env_cmd=None):
        self.queue = queue
        self.nprocs = nprocs
        self.env_cmd = env_cmd

    def refresh_preferred_queue(self) -> bool:
        return True

    def generate_sh(self, workdir, commands: [str], name):
        pass

    def submit(self):
        pass

    def is_running(self, name):
        pass

    def kill_job(self, name):
        pass

    def get_all_jobs(self) -> [Job]:
        return []

    @property
    def n_running_jobs(self) -> int:
        return len(self.get_all_jobs())