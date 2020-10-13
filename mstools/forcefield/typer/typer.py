from ... import logger


class Typer():
    '''
    Base class for typing engines.

    Typing engine assigns atom types for atoms in a molecule or topology by some predefined rules.
    It is the very first step for force field assignment.
    It is also the basis of force field development.
    A well defined typing rule will make the force field development much less painful.
    '''

    def type(self, topology):
        '''
        Assign types for all atoms in the topology.

        The typer will iterate all the molecules belong to the topology and call :func:`type_molecule` on each of them.
        The :attr:`~mstools.topology.Atom.type` attribute of all atoms in the topology will be updated.

        Parameters
        ----------
        topology : Topology
        '''
        for mol in topology.molecules:
            self.type_molecule(mol)

    def type_molecule(self, molecule):
        '''
        Assign types for all the atoms in the molecule.
        This method should be implemented by subclasses.

        Parameters
        ----------
        molecule : Molecule
        '''
        raise NotImplementedError('This method haven\'t been implemented')
