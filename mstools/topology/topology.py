import numpy as np
import copy
from .atom import Atom
from .connectivity import *
from .molecule import Molecule
from .unitcell import UnitCell
from ..forcefield import FFSet


class Topology():
    def __init__(self):
        self.remark = ''
        # TODO Better to set the default cell to None
        self.cell = UnitCell([0, 0, 0])
        self._molecules: [Molecule] = []
        self._atoms: [Atom] = []

    def __deepcopy__(self, memodict={}):
        top = Topology()
        top.remark = self.remark
        top.cell = UnitCell(self.cell.vectors)
        top.init_from_molecules(self._molecules, deepcopy=True)
        return top

    def init_from_molecules(self, molecules: [Molecule], numbers=None, deepcopy=False):
        '''
        Initialize a topology from a bunch of molecule.
        The molecules and atoms are deep copied if numbers is not None or deepcopy is True.
        This method should be called as long as atoms or molecules are added to or removed from topology.
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
        '''
        Assign the id of molecules and atoms belongs to this topology
        '''
        for i, mol in enumerate(self._molecules):
            mol.id = i
        for i, atom in enumerate(self._atoms):
            atom.id = i

    def add_molecule(self, molecule: Molecule):
        '''
        Add molecule into topology.
        Note that the molecule and atoms are not deep copied, the reference are passed.
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
        '''
        positions should always be numpy array
        '''
        return np.array([atom.position for atom in self.atoms])

    @property
    def is_drude(self):
        return any(atom.is_drude for atom in self._atoms)

    @staticmethod
    def open(file, **kwargs):
        from .psf import Psf
        from .pdb import Pdb
        from .lammps import LammpsData
        from .xyz import XyzTopology
        from .zmat import Zmat

        if file.endswith('.psf'):
            return Psf(file, **kwargs)
        if file.endswith('.pdb'):
            return Pdb(file, **kwargs)
        elif file.endswith('.lmp'):
            return LammpsData(file, **kwargs)
        elif file.endswith('.xyz'):
            return XyzTopology(file, **kwargs)
        elif file.endswith('.zmat'):
            return Zmat(file, **kwargs)
        else:
            raise Exception('Unsupported format')

    def write(self, file):
        from .psf import Psf
        from .pdb import Pdb
        from .xyz import XyzTopology

        if file.endswith('.psf'):
            Psf.save_to(self, file)
        elif file.endswith('.pdb'):
            Pdb.save_to(self, file)
        elif file.endswith('.xyz'):
            XyzTopology.save_to(self, file)
        else:
            raise Exception('Unsupported format for topology output')

    @staticmethod
    def save_to(top, file):
        '''
        This method should be implemented by subclasses
        '''
        raise NotImplementedError('This method haven\'t been implemented')

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

    def get_12_13_14_pairs(self) -> ([(Atom, Atom)], [(Atom, Atom)], [(Atom, Atom)]):
        '''
        The pairs only concerns real atoms. Drude particles  will be ignored.
        '''
        pair_12_set = set()
        pair_13_set = set()
        pair_14_set = set()
        for bond in filter(lambda x: not x.is_drude, self.bonds):
            a2, a3 = bond.atom1, bond.atom2
            pair = tuple(sorted([a2, a3]))
            pair_12_set.add(pair)

            for a1 in filter(lambda x: not x.is_drude, a2.bond_partners):
                if a1 != a3:
                    pair = tuple(sorted([a1, a3]))
                    pair_13_set.add(pair)
            for a4 in filter(lambda x: not x.is_drude, a3.bond_partners):
                if a2 != a4:
                    pair = tuple(sorted([a2, a4]))
                    pair_13_set.add(pair)

            for a1 in filter(lambda x: not x.is_drude, a2.bond_partners):
                for a4 in filter(lambda x: not x.is_drude, a3.bond_partners):
                    if a1 != a3 and a2 != a4 and a1 != a4:
                        pair = tuple(sorted([a1, a4]))
                        pair_14_set.add(pair)

        pair_12_list = list(sorted(pair_12_set))
        pair_13_list = list(sorted(pair_13_set - pair_12_set))
        pair_14_list = list(sorted(pair_14_set - pair_13_set.union(pair_12_set)))
        return pair_12_list, pair_13_list, pair_14_list

    def get_drude_pairs(self) -> [(Atom, Atom)]:
        '''
        [(parent, drude)]
        '''
        return [pair for mol in self._molecules for pair in mol.get_drude_pairs()]

    def generate_angle_dihedral_improper(self):
        for mol in self._molecules:
            mol.generate_angle_dihedral_improper()

    def guess_connectivity_from_ff(self, params: FFSet, **kwargs):
        for mol in self._molecules:
            mol.guess_connectivity_from_ff(params, **kwargs)

    def generate_drude_particles(self, params: FFSet, **kwargs):
        for mol in self._molecules:
            mol.generate_drude_particles(params, **kwargs)

    def remove_drude_particles(self):
        for mol in self._molecules:
            mol.remove_drude_particles()

    def assign_mass_from_ff(self, params: FFSet):
        for mol in self._molecules:
            mol.assign_mass_from_ff(params)

    def assign_charge_from_ff(self, params: FFSet):
        for mol in self._molecules:
            mol.assign_charge_from_ff(params)
