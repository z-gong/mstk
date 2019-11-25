import numpy as np
class Atom():
    def __init__(self):
        self.name = ''
        self.type = ''
        self.charge = 0.
        self.position = np.array([0, 0, 0], dtype=float)

class Molecule():
    def __init__(self):
        self.atoms:[Atom] = []
