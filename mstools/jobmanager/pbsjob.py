class PbsJob():
    def __init__(self, id, name, state, workdir):
        self.id = id
        self.name = name
        self.state = state
        self.workdir = workdir
        self.state_str = str(self.state)

    def __repr__(self):
        return '<PbsJob: %i %s %i>' % (self.id, self.name, self.state)

    def __eq__(self, other):
        return self.id == other.id

    class State:
        PENDING = 0
        RUNNING = 1
        DONE = 9

