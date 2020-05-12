from ... import logger


class TypingNotSupportedError(Exception):
    pass


class TypingUndefinedError(Exception):
    pass


class Typer():
    def __init__(self):
        pass

    def type(self, topology):
        for mol in topology.molecules:
            try:
                self.type_molecule(mol)
            except TypingNotSupportedError as e:
                logger.warning('%s not supported by %s: %s' % (
                    str(mol), self.__class__.__name__, str(e)))
            except TypingUndefinedError as e:
                logger.warning('%s not fully typed by %s: %s' % (
                    str(mol), self.__class__.__name__, str(e)))

    def type_molecule(self, molecule):
        '''
        Should be implemented by inheritor
        '''
        raise Exception('This method haven\'t been implemented')
