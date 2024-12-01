import os
import shutil
import numpy as np
import copy
import tempfile
from .atom import Atom
from .molecule import Molecule
from .unitcell import UnitCell
from ..forcefield import ForceField


class Topology():
    '''
    A topology is a combination of molecules and unit cell.

    Parameters
    ----------
    molecules : list of Molecule, optional
        If molecules not provided, the topology should be initialized by calling update_molecules().
    numbers : list of int, optional
        The number of each molecule you want to put in this topology.
    cell : UnitCell, optional
        If cell not provided, a default cell with zero volume will be assumed.

    Notes
    -----
    All the molecules will be deep copied if there are identical molecules in the list, or if numbers is not None.
    The cell will always be deep copied if provided.

    Attributes
    ----------
    remark : str
        A short remark for this topology. Only used for output topology file
    cell : UnitCell
        Unit cell of this topology. Required for output simulation files and guessing connectivity
    '''

    _class_map = {}

    @staticmethod
    def register_format(extension, cls):
        Topology._class_map[extension] = cls

    def __init__(self, molecules=None, numbers=None, cell=None):
        self.remark = ''
        self._molecules: [Molecule] = []
        self._atoms: [Atom] = []
        if molecules is not None:
            self.update_molecules(molecules, numbers)
        if cell is not None:
            self.cell = UnitCell(cell.vectors)
        else:
            self.cell = UnitCell()

    def __repr__(self):
        return f'<Topology: {self.n_molecule} molecules {self.n_residue} residues {self.n_atom} atoms>'

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
        numbers : list of int, Optional
        deepcopy : bool

        Notes
        -----
        The molecules and atoms are deep copied if numbers is not None or deepcopy is True.
        The molecules and atoms are deep copied if there are identical molecules in the list.
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

        # cache atoms
        self._atoms = [atom for mol in self._molecules for atom in mol.atoms]

        # assign the id of molecules, atoms and residues belonging to this topology
        for i, mol in enumerate(self._molecules):
            mol.id = i
        for i, atom in enumerate(self._atoms):
            atom.id = i
        for i, residue in enumerate(self.residues):
            residue.id = i

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
            atom.position = positions[i]

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

    def scale_with_packmol(self, numbers, packmol=None, seed=1, tempdir=None):
        '''
        Enlarge the number of molecules in this topology with Packmol.

        It is used for building real simulation box from several molecules.
        It generates molecules with appropriate positions and adds them into the topology.
        The unit cell should already be defined, and it will not be touched during the scaling.

        An input file for packmol and xyz files for molecules will be created under `tempdir`.
        If packmol is provided, then it will be invoked to pack the simulation box,
        and the generated positions will be loaded back to the topology,
        and the `tempdir` will be removed.
        Otherwise, the topology will be enlarged without correct positions,
        and the files under `tempdir` will be kept so that positions can be built manually and loaded back to the topology later.

        Parameters
        ----------
        numbers : int or list of int
            If numbers is a int, the topology will be scaled by this times.
            If numbers is a list of int, then each molecule will be scaled by corresponding times.
        packmol : Packmol, Optional
            If is None, will generate input files for packmol.
        seed : int
            Seed for random number for Packmol
        tempdir : str, Optional
            The temporary directory to keep intermediate files for packmol. It should exist already
            If is None, will create a directory under system temporary directory
        '''
        if type(numbers) is int:
            numbers = [numbers] * self.n_molecule
        else:
            if len(numbers) != self.n_molecule:
                raise Exception('Length of numbers should equal to molecules')
            if any(type(n) is not int for n in numbers):
                raise Exception('Numbers should be list of integers')

        if self.cell.volume == 0:
            raise Exception('Unit cell required by Packmol')

        if not self.cell.is_rectangular:
            raise Exception('Triclinic unit cell haven\'t been implemented')

        if not self.has_position:
            raise Exception('Positions required by Packmol')

        from ..wrapper.packmol import Packmol
        if packmol is not None and not isinstance(packmol, Packmol):
            raise Exception('Invalid Packmol')

        tempdir = tempdir or tempfile.mkdtemp(prefix='mstk-packmol-')

        xyz_files = []
        for i, mol in enumerate(self._molecules):
            xyz = tempdir + '/_MO_%i.xyz' % i
            xyz_files.append(xyz)
            Topology([mol]).write(xyz)

        tmp_inp = tempdir + '/_pack.inp'
        tmp_out = tempdir + '/_out.xyz'
        box = self.cell.get_size()
        if packmol is not None:
            packmol.build_box(xyz_files, numbers, tmp_out, size=box - 0.2, inp_file=tmp_inp, seed=seed, silent=True)
            try:
                xyz = Topology.open(tmp_out)
            except:
                raise Exception('Building failed. xyz not found')

            if xyz.n_atom != sum(mol.n_atom * n for mol, n in zip(self._molecules, numbers)):
                raise Exception('Building failed. incorrectly information in xyz')

            self.update_molecules(self._molecules, numbers)
            self.set_positions(xyz.positions)
            shutil.rmtree(tempdir)
        else:
            Packmol.gen_inp(xyz_files, numbers, tmp_out, size=box - 0.2, inp_file=tmp_inp, seed=seed)

    def scale_box(self, scale, rigid_group='atom'):
        '''
        Scale the box and positions of atoms

        Parameters
        ----------
        scale : float or np.array of shape (3,)
        rigid_group : ['atom', 'molecule', 'residue']
            If set to atom, the positions of all atoms will be equally scaled
            If set to molecule, each molecule will be treated as a rigid body during scaling
            If set to residue, each residue will be treated as a rigid body during scaling
        '''
        self.cell.set_box(self.cell.vectors * scale)
        if rigid_group == 'atom':
            for atom in self.atoms:
                atom.position = atom.position * scale
                return None

        if rigid_group == 'molecule':
            groups = self.molecules
        elif rigid_group == 'residue':
            groups = self.residues
        else:
            raise Exception(f'Invalid rigid_group: {rigid_group}')

        positions = self.positions
        masses = np.array([a.mass for a in self.atoms])
        for group in groups:
            ids = [a.id for a in group.atoms]
            pos_com = np.sum(positions[ids] * masses[ids, np.newaxis], axis=0) / np.sum(masses[ids])
            for atom in group.atoms:
                atom.position = pos_com * scale + atom.position - pos_com

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
    def n_residue(self):
        '''
        Number of residues belong to this topology

        Returns
        -------
        n : int

        '''
        return sum(mol.n_residue for mol in self._molecules)

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
    def residues(self):
        '''
        List of residues belonging to this topology

        The residue list is generated on the fly because it is called irregularly.
        Therefore, the performance of this method is not of concern.

        Returns
        -------
        residues : list of Residue
        '''
        return [residue for mol in self._molecules for residue in mol.residues]

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

    @property
    def virtual_site_classes(self):
        '''
        A set of all the VirtualSite classes in the topology

        Returns
        -------
        classes : set of subclass of VirtualSite
        '''
        return set(type(atom.virtual_site) for atom in self._atoms if atom.virtual_site is not None)

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

        Returns
        -------
        topology : Topology
        '''
        if file.startswith(':'):
            mol = Molecule.from_smiles(file[1:])
            return Topology([mol])

        basename, ext = os.path.splitext(file)
        ext = ext.lower()
        try:
            cls = Topology._class_map[ext]
        except KeyError:
            raise Exception('Unsupported topology format: ' + file)

        return cls(file, **kwargs).topology

    def write(self, file, **kwargs):
        '''
        Write topology to a file.

        The format will be determined from the extension of the file.

        Parameters
        ----------
        file : str
        '''
        basename, ext = os.path.splitext(file)
        ext = ext.lower()
        try:
            cls = Topology._class_map[ext]
        except KeyError:
            raise Exception('Unsupported format for topology: ' + ext)

        cls.save_to(self, file, **kwargs)

    @staticmethod
    def save_to(top, file, **kwargs):
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
            from openmm import app as omm_app
            from openmm.app.element import Element as omm_Element
        except ImportError:
            raise Exception('cannot import openmm')

        omm_top = omm_app.Topology()
        omm_chain = omm_top.addChain()
        d_omm_residues = {}
        for residue in self.residues:
            omm_res = omm_top.addResidue(residue.name, omm_chain)
            for atom in residue.atoms:
                d_omm_residues[atom] = omm_res
        for atom in self.atoms:
            try:
                omm_element = omm_Element.getBySymbol(atom.symbol)
            except:
                omm_element = None
            omm_top.addAtom(atom.name, omm_element, d_omm_residues[atom])

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
        from openmm.app import Topology as ommTopology, Residue as ommResidue, Atom as ommAtom
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

    def get_virtual_site_pairs(self):
        '''
        Retrieve all the virtual site pairs belong to this topology

        TODO It assumes no more than one virtual site is attached to each atom

        Returns
        -------
        pairs :  list of tuple of Atom
            [(parent, atom_virtual_site)]

        See Also
        --------
        Molecule.get_virtual_site_pairs
        '''
        return [pair for mol in self._molecules for pair in mol.get_virtual_site_pairs()]

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

    def generate_angle_dihedral_improper(self, dihedral=True, improper=True):
        '''
        Generate angles, dihedrals and impropers from bonds.

        See Also
        --------
        Molecule.generate_angle_dihedral_improper
        '''
        for mol in self._molecules:
            mol.generate_angle_dihedral_improper(dihedral=dihedral, improper=improper)

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
        Generate Drude particles from DrudePolarTerms in a polarizable force field

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

    def generate_virtual_sites(self, ff, **kwargs):
        '''
        Generate virtual sites from VirtualSiteTerms in a force field

        Parameters
        ----------
        ff : ForceField
        kwargs : dict

        See Also
        --------
        Molecule.generate_virtual_sites
        '''
        mol: Molecule
        for mol in self._molecules:
            mol.generate_virtual_sites(ff, update_topology=False)
        self.update_molecules(self._molecules, deepcopy=False)

    def remove_virtual_sites(self):
        '''
        Remove all the virtual sites from this topology

        It is useful for build non-virtual-site model from virtual-site model

        See Also
        --------
        Molecule.remove_virtual_sites
        '''
        mol: Molecule
        for mol in self._molecules:
            mol.remove_virtual_sites(update_topology=False)
        self.update_molecules(self._molecules, deepcopy=False)
