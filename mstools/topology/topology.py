import numpy as np
import copy
from io import IOBase
from .molecule import Atom, Bond, Angle, Dihedral, Improper, Molecule
from .unitcell import UnitCell


class Topology():
    def __init__(self):
        self.remark = ''
        self.cell = UnitCell([0, 0, 0])
        self._molecules: [Molecule] = []
        self._atoms: [Atom] = []
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
        self.cell = UnitCell(topology.cell.vectors)
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
    def positions(self):
        return [atom.position for atom in self.atoms]

    @property
    def is_drude(self):
        return any(atom.is_drude for atom in self._atoms)

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

    def to_omm_topology(self):
        try:
            from simtk.openmm import app as omm_app
            from simtk.openmm.app.element import Element as omm_Element
        except ImportError:
            raise Exception('cannot import openmm')

        omm_top = omm_app.Topology()
        omm_chain = omm_top.addChain()
        for mol in self._molecules:
            omm_residue = omm_top.addResidue(mol.name, omm_chain)
            for atom in mol.atoms:
                try:
                    omm_element = omm_Element.getBySymbol(atom.symbol)
                except:
                    omm_element = None
                omm_top.addAtom(atom.name, omm_element, omm_residue)

        omm_top.setPeriodicBoxVectors(self.cell.vectors)
        omm_atoms = list(omm_top.atoms())
        for bond in self.bonds:
            omm_top.addBond(omm_atoms[bond.atom1.id], omm_atoms[bond.atom2.id])

        return omm_top

    def get_12_13_14_exclusions(self) -> ([Atom], [Atom], [Atom]):
        pair_12_set = set()
        pair_13_set = set()
        pair_14_set = set()
        for bond in self.bonds:
            a2, a3 = bond.atom1, bond.atom2
            pair = tuple(sorted([a2, a3]))
            pair_12_set.add(pair)

            for a1 in a2.bond_partners:
                if a1 != a3:
                    pair = tuple(sorted([a1, a3]))
                    pair_13_set.add(pair)
            for a4 in a3.bond_partners:
                if a2 != a4:
                    pair = tuple(sorted([a2, a4]))
                    pair_13_set.add(pair)

            for a1 in a2.bond_partners:
                for a4 in a3.bond_partners:
                    if a1 != a3 and a2 != a4 and a1 != a4:
                        pair = tuple(sorted([a1, a4]))
                        pair_14_set.add(pair)

        pair_12_list = list(sorted(pair_12_set))
        pair_13_list = list(sorted(pair_13_set - pair_12_set))
        pair_14_list = list(sorted(pair_14_set - pair_13_set.union(pair_12_set)))
        return pair_12_list, pair_13_list, pair_14_list
