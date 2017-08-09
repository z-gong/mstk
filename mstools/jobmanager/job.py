class JobState():
    PENDING = 0
    RUNNING = 1
    DONE = 9


class Job():
    def __init__(self, id, name, state, workdir):
        self.id = id
        self.name = name
        self.state = state
        self.workdir = workdir

    def __repr__(self):
        return '<PBS Job: %i %s %i>' % (self.id, self.name, self.state)
