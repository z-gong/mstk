import datetime
import os
import pwd
import time
from .pbsjob import PbsJob


class JobManager:
    def __init__(self, queue=None, nprocs=1, env_cmd=None):
        self.queue = queue
        self.nprocs = nprocs
        self.env_cmd = env_cmd or ''
        self.username = pwd.getpwuid(os.getuid()).pw_name

        self.stored_jobs_expire = 60  # seconds
        self.stored_jobs = []
        self.updated = False
        self.last_update = None

    @property
    def all_jobs(self) -> [PbsJob]:
        if not self.updated or (datetime.datetime.now() - self.last_update).total_seconds() >= self.stored_jobs_expire:
            self.update_stored_jobs()
        return self.stored_jobs

    def update_stored_jobs(self):
        jobs = []
        for i in range(3, 0, -1):
            # in case get_all_jobs() raise Exception
            try:
                jobs += self.get_all_jobs()
            except Exception as e:
                print(repr(e))
            if i > 1:
                time.sleep(1)
        for job in reversed(jobs):  # reverse the job list, in case jobs with same name
            if job not in self.stored_jobs:
                self.stored_jobs.append(job)

        self.stored_jobs.reverse()  # reverse the job list
        self.last_update = datetime.datetime.now()

    def refresh_preferred_queue(self) -> bool:
        return True

    def generate_sh(self, workdir, commands: [str], name):
        pass

    def submit(self) -> bool:
        pass

    def get_id_from_name(self, name: str) -> int:
        for job in self.all_jobs:
            if job.name == name:
                return job.id
        return None

    def is_running(self, name) -> bool:
        for job in self.all_jobs:
            if job.name == name:
                if job.state in [PbsJob.State.PENDING, PbsJob.State.RUNNING]:
                    return True
                else:
                    return False
        return False

    def kill_job(self, name):
        pass

    def get_all_jobs(self) -> [PbsJob]:
        return []

    @property
    def n_running_jobs(self) -> int:
        self.update_stored_jobs()  # must update stored jobs first
        n = 0
        for job in self.all_jobs:
            if job.state != PbsJob.State.DONE:
                n += 1
        return n
