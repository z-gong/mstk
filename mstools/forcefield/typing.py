class TypingEngine():
    def __init__(self):
        pass

    def type(self, topology):
        for mol in topology.molecules:
            self.type_molecule(mol)

    def type_molecule(self, molecule):
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
