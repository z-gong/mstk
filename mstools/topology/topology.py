import shutil
import numpy as np
import copy
import tempfile
from .atom import Atom
from .connectivity import *
from .molecule import Molecule
from .unitcell import UnitCell
from ..forcefield import ForceField
from .. import logger


class Topology():
    '''
    A topology is a combination of molecules and unit cell.
    '''

    def __init__(self, molecules=None, cell=None):
        '''
        Initialize a topology from molecules and unit cell.

        If molecules not provided, the topology should be initialized by calling update_molecules().
        If cell not provided, a default cell with zero volume will be assumed.

        Parameters
        ----------
        molecules : list of Molecule
        cell : UnitCell

        Notes
        -----
        If there are identical molecules provided, then all the molecules will be deep copied.
        The cell will always be deep copied.
        '''

        #: remark for this topology. Only used for output topology file
        self.remark = ''
        #: unit cell of this topology. Required for output simulation files and guessing connectivity
        self.cell = UnitCell([0, 0, 0])
        self._molecules: [Molecule] = []
        self._atoms: [Atom] = []
        if molecules is not None:
            self.update_molecules(molecules)
        if cell is not None:
            self.cell.set_box(cell.vectors)

    def __deepcopy__(self, memodict={}):
        top = Topology()
        top.remark = self.remark
        top.cell = UnitCell(self.cell.vectors)
        top.update_molecules(self._molecules, deepcopy=True)
        return top

    def update_molecules(self, molecules, numbers=None, deepcopy=False):
        '''
        Update the molecules contained in this topology.

        The index of molecules and atoms will be updated.
        Therefore, this method should be called as long as atoms or molecules are added to or removed from topology.

        Parameters
        ----------
        molecules : list of Molecule
        numbers : list of int
        deepcopy : bool

        Notes
        -----
        The molecules and atoms are deep copied if numbers is not None or deepcopy is True.
        The molecules and atoms are deep copied if there are identical molecules provided.
        '''
        if not isinstance(molecules, (list, tuple)) \
                or any(type(mol) is not Molecule for mol in molecules):
            raise Exception('A list of molecules should be provided')

        if len(set(molecules)) < len(molecules) and numbers is None:
            numbers = [1] * len(molecules)

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

        self._assign_id()

    def _assign_id(self):
        '''
        Assign the id of molecules and atoms belongs to this topology
        '''
        for i, mol in enumerate(self._molecules):
            mol.id = i
        for i, atom in enumerate(self._atoms):
            atom.id = i

    def add_molecule(self, molecule):
        '''
        Add a molecule into this topology.

        The index of the added molecule and atoms will be updated.

        Parameters
        ----------
        molecule : Molecule

        Notes
        -----
        The molecule and atoms are not deep copied, the references are passed.

        '''
        molecule.id = self.n_molecule
        molecule._topology = self
        self._molecules.append(molecule)
        for atom in molecule.atoms:
            atom.id = self.n_atom
            self._atoms.append(atom)

    def set_positions(self, positions):
        '''
        Set the positions of all atoms in this topology

        Parameters
        ----------
        positions : array_like
        '''
        if self.n_atom != len(positions):
            raise Exception('Length of positions should equal to the number of atoms')
        for i, atom in enumerate(self.atoms):
            atom.position = np.array(positions[i])

    def get_unique_molecules(self, deepcopy=True):
        '''
        Get the unique molecules and the number of molecules that are similar to these unique molecules.

        This is mainly used for exporting GROMACS topol file.
        By default, the unique molecules are deep copied.

        Parameters
        ----------
        deepcopy : bool

        Returns
        -------
        molecule_number : dict [Molecule, int]
        '''
        mols_unique = []
        flags = [-1] * self.n_molecule
        mol_last: Molecule = None
        for i, mol in enumerate(self._molecules):
            if mol_last is not None and mol.is_similar_to(mol_last):
                flags[mol.id] = flags[mol_last.id]
            else:
                mols_unique.append(mol)
                flags[mol.id] = mol.id
            mol_last = mol

        if deepcopy:
            return {copy.deepcopy(mol): flags.count(mol.id) for mol in mols_unique}
        else:
            return {mol: flags.count(mol.id) for mol in mols_unique}

    def scale_with_packmol(self, numbers, packmol=None):
        '''
        Enlarge the number of molecules in this topology with Packmol.

        It is used for building real simulation box from several different molecules.
        It generates molecules with appropriate positions and addes them into the topology.
        The unit cell should already be defined, and it will not be touched during the scaling.

        If packmol is provided, then it will be invoked to pack the simulation box,
        and the generated positions will be loaded back to the topology.
        Otherwise, the topology will be enlarged without correct positions,
        and a input file for packmol and xyz files for molecules will be saved in current folder
        so that positions can be built manually and loaded back to the topology later

        Parameters
        ----------
        numbers : int or list of int
            If numbers is a int, the topology will be scaled by this times.
            If numbers is a list of int, the each molecule will be scaled by corresponding times.
        packmol : Packmol
        '''
        if type(numbers) is int:
            numbers = [numbers] * self.n_molecule
        else:
            if len(numbers) != self.n_molecule:
                raise Exception('Length of numbers should equal to molecules')
            if any(type(n) is not int for n in numbers):
                raise Exception('Numbers should be list of integers')

        if self.cell.volume == 0:
            raise Exception('Unit cell should be defined first')

        if not self.cell.is_rectangular:
            raise Exception('Triclinic unit cell haven\'t been implemented')

        if not self.has_position:
            raise Exception('Positions are required for scaling box')

        from ..wrapper.packmol import Packmol
        packmol: Packmol

        xyz_files = []
        for i, mol in enumerate(self._molecules):
            _top = Topology([mol])
            if packmol is not None:
                xyz = tempfile.NamedTemporaryFile(suffix='.xyz', prefix='in-').name
            else:
                xyz = '_MO_%i.xyz' % i
            _top.write(xyz)
            xyz_files.append(xyz)

        if packmol is not None:
            print('Build with Packmol ...')
            tmp_inp = tempfile.NamedTemporaryFile(suffix='.inp', prefix='pack-').name
            tmp_out = tempfile.NamedTemporaryFile(suffix='.xyz', prefix='out-').name
            packmol.build_box(xyz_files, numbers, size=self.cell.size * 10 - 2.0, output=tmp_out, inp_file=tmp_inp,
                              silent=True)
            try:
                xyz = Topology.open(tmp_out)
            except:
                raise Exception('Building failed. xyz not found')

            if xyz.n_atom != sum(mol.n_atom * numbers[i] for i, mol in enumerate(self._molecules)) \
                    or not xyz.has_position:
                raise Exception('Building failed. incorrectly information in xyz')

            self.update_molecules(self._molecules, numbers)
            self.set_positions(xyz.positions)
        else:
            tmp_inp = '_pack.inp'
            tmp_out = '_out.xyz'
            Packmol.gen_inp(xyz_files, numbers, size=self.cell.size * 10 - 2.0, output=tmp_out, inp_file=tmp_inp)

    @property
    def n_molecule(self):
        '''
        Number of molecules belong to this topology

        Returns
        -------
        n : int
        '''
        return len(self._molecules)

    @property
    def n_atom(self):
        '''
        Number of atoms belong to this topology

        Returns
        -------
        n : int
        '''
        return len(self._atoms)

    @property
    def molecules(self):
        '''
        List of molecules belong to this topology

        Returns
        -------
        molecules : list of Molecule
        '''
        return self._molecules

    @property
    def atoms(self):
        '''
        List of atoms belong to this topology

        Returns
        -------
        atoms: list of Atom
        '''
        return self._atoms

    @property
    def n_bond(self):
        '''
        Number of bonds belong to this topology

        Returns
        -------
        n : int
        '''
        return sum([mol.n_bond for mol in self._molecules])

    @property
    def n_angle(self):
        '''
        Number of angles belong to this topology

        Returns
        -------
        n : int
        '''
        return sum([mol.n_angle for mol in self._molecules])

    @property
    def n_dihedral(self):
        '''
        Number of dihedrals belong to this topology

        Returns
        -------
        n : int
        '''
        return sum([mol.n_dihedral for mol in self._molecules])

    @property
    def n_improper(self):
        '''
        Number of impropers belong to this topology

        Returns
        -------
        n : int
        '''
        return sum([mol.n_improper for mol in self._molecules])

    @property
    def bonds(self):
        '''
        List of bonds belong to this topology

        Returns
        -------
        bonds : list of Bond
        '''
        return [bond for mol in self._molecules for bond in mol.bonds]

    @property
    def angles(self):
        '''
        List of angles belong to this topology

        Returns
        -------
        angles : list of Angle
        '''
        return [angle for mol in self._molecules for angle in mol.angles]

    @property
    def dihedrals(self):
        '''
        List of dihedrals belong to this topology

        Returns
        -------
        dihedrals : list of Dihedral
        '''
        return [dihedral for mol in self._molecules for dihedral in mol.dihedrals]

    @property
    def impropers(self):
        '''
        List of impropers belong to this topology

        Returns
        -------
        impropers : list of Improper
        '''
        return [improper for mol in self._molecules for improper in mol.impropers]

    @property
    def has_position(self):
        '''
        Whether or not all the atoms in the topology have positions

        Returns
        -------
        has : bool
        '''
        return all(atom.has_position for atom in self.atoms)

    @property
    def positions(self):
        '''
        Positions of all the atoms in this topology

        Returns
        -------
        positions : array_like
        '''
        return np.array([atom.position for atom in self.atoms])

    @property
    def is_drude(self):
        '''
        Whether or not this topology represents a Drude model

        Returns
        -------
        is : bool
        '''
        return any(atom.is_drude for atom in self._atoms)

    @property
    def has_virtual_site(self):
        '''
        Whether or not this topology contains virtual site

        Returns
        -------
        has : bool
        '''
        return any(atom.virtual_site is not None for atom in self._atoms)

    @staticmethod
    def open(file, **kwargs):
        '''
        Load topology from a file or smiles string

        The format will be determined from the extension of the file.

        Parameters
        ----------
        file : str
            If start with ':', then it will be treated as SMLIES string.
            Otherwise, it is the file to be read.
        kwargs : dict, optional

        Returns
        -------
        topology : Topology
        '''
        from .psf import Psf
        from .pdb import Pdb
        from .lammps import LammpsData
        from .xyz import XyzTopology
        from .zmat import Zmat

        if file.startswith(':'):
            mol = Molecule.from_smiles(file[1:])
            return Topology([mol])
        elif file.endswith('.psf'):
            return Psf.read(file, **kwargs)
        elif file.endswith('.pdb'):
            return Pdb.read(file, **kwargs)
        elif file.endswith('.lmp'):
            return LammpsData.read(file, **kwargs)
        elif file.endswith('.xyz'):
            return XyzTopology.read(file, **kwargs)
        elif file.endswith('.zmat'):
            return Zmat.read(file, **kwargs)
        else:
            raise Exception('Unsupported format')

    def write(self, file):
        '''
        Write topology to a file.

        The format will be determined from the extension of the file.

        Parameters
        ----------
        file : str
        '''
        from . import Psf, Pdb, XyzTopology

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
        '''
        Generate OpenMM topology from this topology

        Returns
        -------
        omm_topology : mm.app.Topology
        '''
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

        if self.cell.volume != 0:
            omm_top.setPeriodicBoxVectors(self.cell.vectors)
        omm_atoms = list(omm_top.atoms())
        for bond in self.bonds:
            omm_top.addBond(omm_atoms[bond.atom1.id], omm_atoms[bond.atom2.id])

        return omm_top

    def init_from_omm_topology(self, omm_top):
        '''
        Init the topology from OpenMM topology.

        This is only used for write trajectory, therefore only the atomic symbol is loaded.
        It should not be used for other purpose.

        Parameters
        ----------
        omm_top : mm.app.Topology
        '''
        from simtk.openmm.app import Topology as ommTopology, Residue as ommResidue, Atom as ommAtom
        molecules = []
        omm_top: ommTopology
        omm_res: ommResidue
        omm_atom: ommAtom
        for omm_res in omm_top.residues():
            mol = Molecule(omm_res.name)
            for omm_atom in omm_res.atoms():
                atom = Atom(omm_atom.name)
                try:
                    atom.symbol = omm_atom.element.symbol
                except:
                    pass
                mol.add_atom(atom)
            molecules.append(mol)
        self.update_molecules(molecules)

    def get_12_13_14_pairs(self):
        '''
        Retrieve all the 1-2, 1-3 and 1-4 pairs based on the bond information.

        The pairs only concerns real atoms. Drude particles  will be ignored.

        Returns
        -------
        pairs12 : list of tuple of Atom
        pairs13 : list of tuple of Atom
        pairs14 : list of tuple of Atom

        See Also
        --------
        Molecule.get_12_13_14_pairs
        '''
        pair_12_list = []
        pair_13_list = []
        pair_14_list = []
        mol: Molecule
        for mol in self._molecules:
            pairs_12, pairs_13, pairs_14 = mol.get_12_13_14_pairs()
            pair_12_list += pairs_12
            pair_13_list += pairs_13
            pair_14_list += pairs_14
        return pair_12_list, pair_13_list, pair_14_list

    def get_drude_pairs(self):
        '''
        Retrieve all the Drude dipole pairs belong to this topology

        Returns
        -------
        pairs :  list of tuple of Atom
            [(parent, drude)]

        See Also
        --------
        Molecule.get_drude_pairs
        '''
        return [pair for mol in self._molecules for pair in mol.get_drude_pairs()]

    def generate_angle_dihedral_improper(self):
        '''
        Generate angles, dihedrals and impropers from bonds.

        See Also
        --------
        Molecule.generate_angle_dihedral_improper
        '''
        for mol in self._molecules:
            mol.generate_angle_dihedral_improper()

    def guess_connectivity_from_ff(self, ff, bond_limit=0.25, bond_tolerance=0.025, angle_tolerance=None, pbc=''):
        '''
        Guess the connectivity between atoms from force field parameters.

        Parameters
        ----------
        ff : ForceField
        bond_limit : float
        bond_tolerance : float
        angle_tolerance : float
        pbc : str

        See Also
        --------
        Molecule.guess_connectivity_from_ff
        '''
        mol: Molecule
        for mol in self._molecules:
            mol.guess_connectivity_from_ff(ff, bond_limit, bond_tolerance,
                                           angle_tolerance, pbc, self.cell)

    def generate_drude_particles(self, ff, **kwargs):
        '''
        Generate Drude particles from DrudeTerms in a polarizable force field

        Parameters
        ----------
        ff : ForceField
        kwargs : dict

        See Also
        --------
        Molecule.generate_drude_particles
        '''
        mol: Molecule
        for mol in self._molecules:
            mol.generate_drude_particles(ff, update_topology=False)
        self.update_molecules(self._molecules, deepcopy=False)

    def remove_drude_particles(self):
        '''
        Remove all the Drude particles from this topology

        It is useful for build non-polarizable model from polarizable model

        See Also
        --------
        Molecule.remove_drude_particles
        '''
        mol: Molecule
        for mol in self._molecules:
            mol.remove_drude_particles(update_topology=False)
        self.update_molecules(self._molecules, deepcopy=False)

    def assign_mass_from_ff(self, ff):
        '''
        Assign masses for all atoms and Drude particles from the AtomType saved in force field.

        Parameters
        ----------
        ff : ForceField

        See Also
        --------
        Molecule.assign_mass_from_ff
        '''
        for mol in self._molecules:
            mol.assign_mass_from_ff(ff)

    def assign_charge_from_ff(self, ff):
        '''
        Assign charges for all atoms and Drude particles from the AtomType saved in force field.

        Parameters
        ----------
        ff : ForceField

        See Also
        --------
        Molecule.assign_charge_from_ff
        '''
        for mol in self._molecules:
            mol.assign_charge_from_ff(ff)
