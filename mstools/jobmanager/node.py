class Node():
    def __init__(self, name, queue, nprocs, n_free_cores, state):
        self.name = name
        self.queue = queue
        self.nprocs = nprocs
        self.n_free_cores = n_free_cores
        self.state = state
