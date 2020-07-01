from ... import logger


class Typer():
    '''
    Typing engine assigns atom types for atoms in a molecule or topology by some predefined rules.
    It is the very first step for force field assignment.
    It is also the basis of force field development.
    A well defined typing rule will make the force field development much less painful.
    '''
    def __init__(self):
        '''
        This method should be implemented by subclasses
        '''
        pass

    def type(self, topology):
        '''
        Type a topology.

        Parameters
        ----------
        topology : Topology

        Returns
        -------

        '''
        for mol in topology.molecules:
            self.type_molecule(mol)

    def type_molecule(self, molecule):
        '''
        Type a molecule.
        This method should be implemented by inheritor

        Parameters
        ----------
        molecule : Molecule

        Returns
        -------
        '''
        raise NotImplementedError('This method haven\'t been implemented')
