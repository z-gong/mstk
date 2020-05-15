from ... import logger


class Typer():
    def __init__(self):
        pass

    def type(self, topology):
        for mol in topology.molecules:
            self.type_molecule(mol)

    def type_molecule(self, molecule):
        '''
        Should be implemented by inheritor
        '''
        raise NotImplementedError('This method haven\'t been implemented')
