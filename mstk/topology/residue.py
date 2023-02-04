class Residue:
    '''
    A residue is a group of consecutive atoms, e.g. the amino acid in protein, the repeating unit in polymer.

    It can be considered as a special connectivity without energy contribution.
    It is mainly for the convenience of visualization of polymers.
    Small molecule usually contains one single residue.
    '''

    def __init__(self, name='UNK'):
        self.id = -1
        self.id_in_mol = -1
        self.name = name
        self._atoms = []

    def __repr__(self):
        return f'<Residue: {self.name} {self.id}>'

    @property
    def n_atom(self):
        return len(self._atoms)

    @property
    def atoms(self):
        return self._atoms

    def _add_atom(self, atom):
        atom._residue = self
        self._atoms.append(atom)

    def _remove_atom(self, atom):
        atom._residue = None
        self._atoms.remove(atom)
