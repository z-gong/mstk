class JobManager:
    def __init__(self, queue=None, nprocs=1):
        self.queue = queue
        self.nprocs = nprocs

    def refresh_preferred_queue(self) -> bool:
        return True

    def generate_sh(self, workdir, commands: [str], name):
        pass

    def submit(self):
        pass

    def get_info(self, name):
        pass

    def kill_job(self, name):
        pass
