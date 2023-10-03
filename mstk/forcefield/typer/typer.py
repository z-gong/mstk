class Typer:
    '''
    Base class for typing engines.

    Typing engine assigns atom types for atoms in a molecule or topology by some predefined rules.
    It is the very first step for force field assignment.
    It is also the basis of force field development.
    A well defined typing rule will make the force field development much less painful.
    '''

    def type(self, top_or_mol):
        '''
        Assign types for all atoms in a topology or molecule.

        The :attr:`~mstk.topology.Atom.type` attribute of all atoms in the topology/molecule will be updated.

        Parameters
        ----------
        top_or_mol: Topology, Molecule
        '''
        from mstk.topology import Topology, Molecule

        if type(top_or_mol) is Topology:
            for mol in top_or_mol.molecules:
                self._type_molecule(mol)
        elif type(top_or_mol) is Molecule:
            self._type_molecule(top_or_mol)
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
