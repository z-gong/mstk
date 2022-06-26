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
        Assign types for all atoms in the topology or molecule.

        The :attr:`~mstk.topology.Atom.type` attribute of all atoms in the topology/molecule will be updated.

        Parameters
        ----------
        topology : Topology, Molecule
        '''
        from mstk.topology import Topology, Molecule

        if type(topology) is Topology:
            for mol in topology.molecules:
                self._type_molecule(mol)
        elif type(topology) is Molecule:
            self._type_molecule(topology)
        else:
            raise Exception('A topology or molecule is expected')

    def _type_molecule(self, molecule):
        '''
        Assign types for all the atoms in the molecule.
        This method should be implemented by subclasses.

        Parameters
        ----------
        molecule : Molecule
        '''
        raise NotImplementedError('This method haven\'t been implemented')
