import numpy as np
import copy
from io import IOBase
from .molecule import Atom, Bond, Angle, Dihedral, Improper, Molecule

class Topology():
    def __init__(self):
        self.remark = ''
        self.is_drude = False
        self._molecules: [Molecule] = []
        self._atoms: [Atom] = []
        self._box = np.array([0, 0, 0], dtype=float)
        self._file = IOBase()

    def __del__(self):
        self.close()

    def close(self):
        self._file.close()

    def init_from_topology(self, topology, deepcopy=False):
        """
        init a new topology from another topology
        this is useful when you want to convert topology into different formats
        by default the molecules are passed as reference without deepcopy
        be careful if you still need to manipulate the original topology

        @param topology: Topology
        @param deepcopy: bool
        """
        self.remark = topology.remark
        self.is_drude = topology.is_drude
        self._box = topology._box
        self.init_from_molecules(self._molecules, deepcopy=deepcopy)

    def init_from_molecules(self, molecules: [Molecule], numbers=None, deepcopy=False):
        '''
        initialize a topology from a bunch of molecules
        the molecules and atoms are deep copied if numbers is not None or deepcopy is True
        '''
        if len(set(molecules)) < len(molecules):
            raise Exception('There are duplicated molecules, consider set the number of molecule')
        if numbers is None and not deepcopy:
            self._molecules = molecules[:]
        else:
            if numbers is None:
                numbers = [1] * len(molecules)
            elif len(molecules) != len(numbers):
                raise Exception('Elements in molecules and numbers do not match')
            self._molecules = []
            for mol, number in zip(molecules, numbers):
                for i in range(number):
                    self._molecules.append(copy.deepcopy(mol))

        for mol in self._molecules:
            mol._topology = self
        self._atoms = [atom for mol in self._molecules for atom in mol.atoms]

        self.assign_id()

    def assign_id(self):
        idx_atom = 0
        for i, mol in enumerate(self._molecules):
            mol.id = i
            for j, atom in enumerate(mol.atoms):
                atom.id = idx_atom
                idx_atom += 1

    def add_molecule(self, molecule: Molecule):
        '''
        Add molecule into topology
        Note that the molecule and atoms are not deep copied, the reference are passed
        '''
        molecule.id = self.n_molecule
        molecule._topology = self
        self._molecules.append(molecule)
        for atom in molecule.atoms:
            atom.id = self.n_atom
            self._atoms.append(atom)

    def set_positions(self, positions):
        if self.n_atom != len(positions):
            raise Exception('Length of positions should equal to the number of atoms')
        for i, atom in enumerate(self.atoms):
            atom.position = np.array(positions[i])
            atom.has_position = True

    def set_velocities(self, velocities):
        if self.n_atom != len(velocities):
            raise Exception('Length of velocities should equal to the number of atoms')
        for i, atom in enumerate(self.atoms):
            atom.velocity = np.array(velocities[i])
            atom.has_velocity = True

    @property
    def n_molecule(self):
        return len(self._molecules)

    @property
    def n_atom(self):
        return len(self._atoms)

    @property
    def molecules(self):
        return self._molecules

    @property
    def atoms(self):
        return self._atoms

    @property
    def n_bond(self):
        return sum([mol.n_bond for mol in self._molecules])

    @property
    def n_angle(self):
        return sum([mol.n_angle for mol in self._molecules])

    @property
    def n_dihedral(self):
        return sum([mol.n_dihedral for mol in self._molecules])

    @property
    def n_improper(self):
        return sum([mol.n_improper for mol in self._molecules])

    @property
    def bonds(self):
        return [bond for mol in self._molecules for bond in mol.bonds]

    @property
    def angles(self):
        return [angle for mol in self._molecules for angle in mol.angles]

    @property
    def dihedrals(self):
        return [dihedral for mol in self._molecules for dihedral in mol.dihedrals]

    @property
    def impropers(self):
        return [improper for mol in self._molecules for improper in mol.impropers]

    @property
    def has_position(self):
        return all(atom.has_position for atom in self.atoms)

    @property
    def has_velocity(self):
        return all(atom.has_velocity for atom in self.atoms)

    @property
    def positions(self):
        return [atom.position for atom in self.atoms]

    @property
    def velocities(self):
        return [atom.velocity for atom in self.atoms]

    @property
    def box(self):
        return self._box

    @box.setter
    def box(self, value):
        if not isinstance(value, (list, tuple, np.ndarray)) or len(value) != 3:
            raise ValueError('box should has three elements')
        self._box = np.array(value)

    @staticmethod
    def open(file, mode='r'):
        from .psf import Psf
        from .lammps import LammpsData
        from .xyz import XyzTopology
        from .zmat import Zmat

        if file.endswith('.psf'):
            return Psf(file, mode)
        elif file.endswith('.lmp'):
            return LammpsData(file, mode)
        elif file.endswith('.xyz'):
            return XyzTopology(file, mode)
        elif file.endswith('.zmat'):
            return Zmat(file, mode)
        else:
            raise Exception('filename for topology not understand')
