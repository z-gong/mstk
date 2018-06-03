class PbsJob():
    def __init__(self, id, name, state, workdir=None, user=None):
        self.id = id
        self.name = name
        self.state = state
        self.workdir = workdir
        self.user = user
        self.state_str = None

    def __repr__(self):
        return '<PbsJob: %i %s %i %s>' % (self.id, self.name, self.state, self.user)

    def __eq__(self, other):
        return self.id == other.id

    class State:
        PENDING = 0
        RUNNING = 1
        DONE = 9
