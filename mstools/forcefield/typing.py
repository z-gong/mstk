from ..topology import Topology, Molecule


class TypingEngine():
    def __init__(self):
        pass

    def type(self, top: Topology):
        for mol in top.molecules:
            self.type_molecule(mol)

    def type_molecule(self, mol: Molecule):
        '''
        Should be implemented by inheritor
        :param mol:
        :return:
        '''
        pass


class DffTypingEngine(TypingEngine):
    def __init__(self, file):
        super().__init__()
        self._file = open(file)

    def type_molecule(self, mol):
        # TODO to be implemented
        pass
