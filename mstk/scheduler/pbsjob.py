class PbsJob:
    '''
    A PbsJob object is a job managed by scheduler.

    Parameters
    ----------
    id : int
        Job ID assigned by job scheduler.
    name : str
        Name of this job in job scheduler.
    state: [0, 1, 9]
        Whether this job is pending (0), running (1) or done (9).
        Successfully finished, Killed and failed are all considered as done.
    workdir : str
        The working directory of this job.
    user : str
        Name of the user who submitted this job.
    queue : str
        The queue this job is submitted to.

    Attributes
    ----------
    id : int
        Job ID assigned by job scheduler.
    name : str
        Name of this job in job scheduler.
    state: [0, 1, 9]
        Whether this job is pending (0), running (1) or done (9).
        Successfully finished, Killed and failed are all considered as done.
    workdir : str
        The working directory of this job.
    user : str
        Name of the user who submitted this job.
    queue : str
        The queue this job is submitted to.
    state_str : str
        The string shown by the job scheduler, which represents the state of the job.
    '''

    def __init__(self, id, name, state, workdir=None, user=None, queue=None):
        self.id = id
        self.name = name
        self.state = state
        self.workdir = workdir
        self.user = user
        self.queue = queue
        self.state_str = None

    def __repr__(self):
        return '<PbsJob: %i %s %i %s>' % (self.id, self.name, self.state, self.user)

    def __eq__(self, other):
        return self.id == other.id

    class State:
        PENDING = 0
        RUNNING = 1
        DONE = 9
