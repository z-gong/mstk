import os
from mstk import MSTK_FORCEFIELD_PATH


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

    @staticmethod
    def open(filename):
        '''
        Load a typer from a type definition file.

        The typing engine is determined by the TypingEngine line in the file.

        Parameters
        ----------
        filename : str
            Type definition file.
            If the file does not exist, will search it under directories defined by MSTK_FORCEFIELD_PATH.

        Returns
        -------
        typer : subclass of Typer
        '''
        from .smarts_typer import SmartsTyper
        from .gaff_typer import GaffTyper

        for dir in ['.'] + MSTK_FORCEFIELD_PATH:
            p = os.path.join(dir, filename)
            if os.path.exists(p):
                filepath = p
                break
        else:
            raise Exception(f'Typing file not found: {filename}')

        engine = 'SmartsTyper'
        with open(p) as f:
            for line in f:
                if line.lower().startswith('typingengine'):
                    engine = line.split()[1]
                    break

        if engine.lower() == 'smartstyper':
            return SmartsTyper(filepath)
        elif engine.lower() == 'gafftyper':
            return GaffTyper(filepath)
        else:
            raise Exception(f'Unknown typing engine: {engine}')
