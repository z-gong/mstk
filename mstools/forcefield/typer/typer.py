from ... import logger


class Typer():
    def __init__(self):
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
