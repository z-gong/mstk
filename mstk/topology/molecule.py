import itertools
import copy
import numpy as np
from mstk import logger
from mstk.chem.rdkit import create_mol_from_smiles
from mstk.chem.element import Element
from mstk.forcefield.ffterm import *
from .atom import Atom
from .virtualsite import *
from .connectivity import *
from .unitcell import UnitCell
from .residue import Residue
from .geometry import find_clusters_in_graph


class Molecule():
    '''
    A molecule is defined as atoms and the connectivity between them.

    The term `molecule` is not strictly a chemical molecule.
    Some atoms may not be connected to any other atoms in the same molecule.
    However, there can not be bonds connecting atoms belong to different molecules.
    Drude particles and virtual sites are also considered as atoms.
    All bond, angles, dihedrals and impropers should be defined explicitly.

    Parameters
    ----------
    name : str

    Attributes
    ----------
    id : int
        Index of this molecule in topology. -1 means information haven\'t been updated by topology
    name : str
        Name of the molecule, not necessarily unique
    '''

    def __init__(self, name='UNK'):
        self.id = -1
        self._name = name
        self._topology = None
        self._atoms: [Atom] = []
        self._bonds: [Bond] = []
        self._angles: [Angle] = []
        self._dihedrals: [Dihedral] = []
        self._impropers: [Improper] = []
        self._rdmol = None  # this is for typing based on SMARTS
        self._is_rdmol_valid = False
        self._default_residue = Residue(name)
        self._added_residues = []

    def __repr__(self):
        return f'<Molecule: {self.name} {self.id}>'

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, val):
        self._name = val
        self._default_residue.name = val

    def __deepcopy__(self, memodict={}):
        '''
        If there are virtual sites, they will be constructed with new atoms also
        If there are residue information, they will also be constructed

        topology will not be copied.
        id will not be copied, because it relies on other molecules in the topology
        '''
        mol = Molecule(self.name)
        for atom in self._atoms:
            mol.add_atom(copy.deepcopy(atom))
        for bond in self._bonds:
            idx1 = bond.atom1.id_in_mol
            idx2 = bond.atom2.id_in_mol
            mol.add_bond(mol._atoms[idx1], mol._atoms[idx2], bond.order)
        for angle in self._angles:
            idx1 = angle.atom1.id_in_mol
            idx2 = angle.atom2.id_in_mol
            idx3 = angle.atom3.id_in_mol
            mol.add_angle(mol._atoms[idx1], mol._atoms[idx2], mol._atoms[idx3])
        for dihedral in self._dihedrals:
            idx1 = dihedral.atom1.id_in_mol
            idx2 = dihedral.atom2.id_in_mol
            idx3 = dihedral.atom3.id_in_mol
            idx4 = dihedral.atom4.id_in_mol
            mol.add_dihedral(mol._atoms[idx1], mol._atoms[idx2], mol._atoms[idx3], mol._atoms[idx4])
        for improper in self._impropers:
            idx1 = improper.atom1.id_in_mol
            idx2 = improper.atom2.id_in_mol
            idx3 = improper.atom3.id_in_mol
            idx4 = improper.atom4.id_in_mol
            mol.add_improper(mol._atoms[idx1], mol._atoms[idx2], mol._atoms[idx3], mol._atoms[idx4])

        # add_atom and add_bond will invalidate rdmol
        if self._is_rdmol_valid:
            from rdkit import Chem
            mol._rdmol = Chem.Mol(self._rdmol)
            mol._is_rdmol_valid = True

        for i, atom in enumerate(self._atoms):
            vsite = atom.virtual_site
            if vsite is not None:
                new_atom = mol.atoms[i]
                new_parents = [mol.atoms[p.id_in_mol] for p in vsite.parents]
                new_atom.virtual_site = VirtualSite.create(vsite.__class__.__name__, new_parents, vsite.parameters)

        for residue in self._added_residues:
            atoms = [mol.atoms[atom.id_in_mol] for atom in residue.atoms]
            mol.add_residue(residue.name, atoms, refresh_residues=False)
        mol.refresh_residues()

        return mol

    @staticmethod
    def from_smiles(smiles):
        '''
        Initialize a molecule from SMILES string.

        RDKit is used for parsing SMILES. The Hydrogen atoms will be created.
        The positions of all atoms will also be automatically generated.
        The SMILES string can contain the name of the molecule at the end, e.g. 'CCCC butane'.

        Parameters
        ----------
        smiles : str

        Returns
        -------
        molecule : Molecule
        '''
        words = smiles.strip().split()
        if not words:
            raise Exception('Invalid SMILES string')

        smiles = words[0]
        if len(words) > 1:
            name = words[1]
        else:
            name = None

        rdmol = create_mol_from_smiles(smiles)
        mol = Molecule.from_rdmol(rdmol, name)

        return mol

    @staticmethod
    def from_rdmol(rdmol, name=None):
        '''
        Initialize a molecule from a RDKit Mol object.
        If the RDKit Mol has conformers, the position of the first conformer will be assigned to the atoms

        Parameters
        ----------
        rdmol : rdkit.Chem.Mol
        name : str
            The name of the molecule. If not provided, the formula will be used as the name.

        Returns
        -------
        molecule : Molecule
        '''
        try:
            from rdkit import Chem
            from rdkit.Chem.rdMolDescriptors import CalcMolFormula
        except ImportError:
            raise ImportError('RDKit not found')

        rdmol = Chem.Mol(rdmol)
        # don't set aromaticity, kekulized bonds are easier to manipulate
        Chem.SanitizeMol(rdmol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY)
        mol = Molecule()
        for i, a in enumerate(rdmol.GetAtoms()):
            atom = Atom()
            element = Element(a.GetAtomicNum())
            atom.name = element.symbol + str(i + 1)
            atom.symbol = element.symbol
            atom.mass = element.mass
            atom.formal_charge = a.GetFormalCharge()
            mol.add_atom(atom)
        if rdmol.GetNumConformers() > 0:
            for atom, pos in zip(mol.atoms, rdmol.GetConformer().GetPositions()):
                atom.position = pos / 10  # convert A to nm
        for b in rdmol.GetBonds():
            atom1 = mol.atoms[b.GetBeginAtomIdx()]
            atom2 = mol.atoms[b.GetEndAtomIdx()]
            d_bond_order = {
                Chem.rdchem.BondType.UNSPECIFIED: Bond.Order.UNSPECIFIED,
                Chem.rdchem.BondType.SINGLE     : Bond.Order.SINGLE,
                Chem.rdchem.BondType.DOUBLE     : Bond.Order.DOUBLE,
                Chem.rdchem.BondType.TRIPLE     : Bond.Order.TRIPLE,
            }
            try:
                order = d_bond_order[b.GetBondType()]
            except KeyError:
                logger.warning('Only single/double/triple/aromatic bond supported. Will discard bond order')
                order = Bond.Order.UNSPECIFIED
            mol.add_bond(atom1, atom2, order)
        mol.generate_angle_dihedral_improper()

        # set aromaticity so that SMARTS matching and is_aromatic works correctly
        Chem.SetAromaticity(rdmol)
        mol._rdmol = rdmol
        mol._is_rdmol_valid = True

        if name is not None:
            mol.name = name
        else:
            mol.name = CalcMolFormula(rdmol)

        return mol

    @property
    def rdmol(self):
        '''
        The `rdkit.Chem.Mol` object associated with this molecule.

        It is required by SmartsTyper typing engine, which performs SMARTS matching on the molecule.
        The `rdmol` attribute will be assigned if the molecule is initialized from SMILES or RDKit Molecule.
        If it is not available, a RDKit molecule will be constructed from atoms and bonds.
        The positions will not be preserved.

        Returns
        -------
        rdmol : rdkit.Chem.Mol
        '''
        try:
            from rdkit import Chem
        except ImportError:
            raise ImportError('RDKit not found')

        if self._is_rdmol_valid:
            return self._rdmol

        if any(b.order == Bond.Order.UNSPECIFIED for b in self.bonds):
            logger.warning(f'Not all bond orders are specified in {self}')

        rwmol = Chem.RWMol()
        for atom in self.atoms:
            rdatom = Chem.Atom(atom.symbol)
            rdatom.SetFormalCharge(atom.formal_charge)
            rdatom.SetNoImplicit(True)  # disable implicit Hs. Otherwise cannot handle radicals
            rwmol.AddAtom(rdatom)
        for bond in self.bonds:
            d_bond_order = {
                Bond.Order.UNSPECIFIED: Chem.rdchem.BondType.UNSPECIFIED,
                Bond.Order.SINGLE     : Chem.rdchem.BondType.SINGLE,
                Bond.Order.DOUBLE     : Chem.rdchem.BondType.DOUBLE,
                Bond.Order.TRIPLE     : Chem.rdchem.BondType.TRIPLE,
            }
            rwmol.AddBond(bond.atom1.id_in_mol, bond.atom2.id_in_mol, d_bond_order[bond.order])
        # set aromaticity so that SMARTS matching and is_aromatic works correctly
        Chem.SanitizeMol(rwmol)
        self._rdmol = rwmol.GetMol()
        self._is_rdmol_valid = True

        return self._rdmol

    def generate_conformers(self, n_conformer=1):
        '''
        Generate several conformers with RDKit.

        The positions will be generated from only elements and bonds. The chiral center will not be respected.

        Parameters
        ----------
        n_conformer : int
            How many conformers to generate

        Returns
        -------
        positions_list : list of array_like
        '''
        try:
            from rdkit.Chem import AllChem as Chem
        except ImportError:
            raise ImportError('RDKit not found')

        positions_list = []
        rdmol = Chem.Mol(self.rdmol)
        useRandomCoords = False
        while len(positions_list) < n_conformer:
            Chem.EmbedMultipleConfs(rdmol, numConfs=n_conformer - len(positions_list), clearConfs=True,
                                    useRandomCoords=useRandomCoords)
            for i in range(rdmol.GetNumConformers()):
                positions_list.append(rdmol.GetConformer(i).GetPositions() / 10)  # A -> nm
            useRandomCoords = True

        return positions_list

    @property
    def topology(self):
        '''
        The topology this molecule belongs to

        Returns
        -------
        topology : Topology

        '''
        return self._topology

    def add_atom(self, atom, residue=None, index=None, update_topology=True):
        '''
        Add an atom to this molecule.

        The id_in_mol attribute of all atoms will be updated after insertion.

        TODO Make residue assignment more robust

        Parameters
        ----------
        atom : Atom
        residue : Residue, Optional
            Add the atom to this residue.
            Make sure the residue belongs to this molecule. For performance concern, this is not checked.
            If set to None, the atom will be added to the default residue.
        index : int, Optional
            If None, the new atom will be the last atom. Otherwise, it will be inserted in front of index-th atom.
        update_topology : bool
            If True, the topology this molecule belongs to will update its atom list and assign id for all atoms and residues.
            Otherwise, you have to re-init the topology manually so that the topological information is correct.
        '''
        atom._molecule = self
        if index is None:
            self._atoms.append(atom)
            atom.id_in_mol = len(self._atoms) - 1
        else:
            self._atoms.insert(index, atom)
            for i in range(index, len(self._atoms)):
                self._atoms[i].id_in_mol = i
        self._is_rdmol_valid = False

        residue = residue or self._default_residue
        residue._add_atom(atom)
        # re-index residues because the residue starts to count
        if residue.n_atom == 1:
            self.refresh_residues(update_topology)

        if self._topology is not None and update_topology:
            self._topology.update_molecules(self._topology.molecules, deepcopy=False)

    def remove_atom(self, atom, update_topology=True):
        '''
        Remove an atom and all the bonds connected to the atom from this molecule.

        The atom will also be removed from its residue.
        The id_in_mol attribute of all atoms will be updated after removal.
        The angle, dihedral and improper involving this atom are untouched.
        Therefore, you may call `generate_angle_dihedral_improper` to refresh the connectivity.

        Parameters
        ----------
        atom : Atom
        update_topology : bool
            If update_topology is True, the topology this molecule belongs to will update its atom list and assign id for all atoms and residues.
            Otherwise, you have to re-init the topology manually so that the topological information is correct.
        '''
        self.remove_atoms([atom], update_topology=update_topology)

    def remove_atoms(self, atoms, update_topology=True):
        '''
        Remove multiple atoms and all the bonds connected to these atoms from this molecule.

        The atom will also be removed from its residue.
        The id_in_mol attribute of all atoms will be updated after removal.
        The angle, dihedral and improper involving this atom are untouched.
        Therefore, you may call `generate_angle_dihedral_improper` to refresh the connectivity.

        Parameters
        ----------
        atom : Atom
        update_topology : bool
            If update_topology is True, the topology this molecule belongs to will update its atom list and assign id for all atoms and residues.
            Otherwise, you have to re-init the topology manually so that the topological information is correct.
        '''
        _refresh_residues = False
        for atom in atoms:
            for bond in atom._bonds[:]:
                self.remove_connectivity(bond)
            # self._atoms.remove(atom)
            atom._molecule = None

            residue = atom.residue
            residue._remove_atom(atom)
            if not _refresh_residues and residue.n_atom == 0:
                _refresh_residues = True

        # it's faster to reconstruct atom list than iteratively calling _atoms.remove()
        self._atoms = list(sorted(set(self._atoms) - set(atoms), key=lambda x: x.id_in_mol))

        for i, at in enumerate(self._atoms):
            at.id_in_mol = i
        self._is_rdmol_valid = False

        # re-index residues because some residues become empty
        if _refresh_residues:
            self.refresh_residues(update_topology)

        if self._topology is not None and update_topology:
            self._topology.update_molecules(self._topology.molecules, deepcopy=False)

    def remove_non_polar_hydrogens(self, update_topology=True):
        '''
        Remove single-coordinated hydrogen atoms bonded to C and Si atoms

        Parameters
        ----------
        update_topology : bool
            If update_topology is True, the topology this molecule belongs to will update its atom list and assign id for all atoms and residues.
            Otherwise, you have to re-init the topology manually so that the topological information is correct.

        Returns
        -------
        ids_removed : list of int
            The list of `id_in_mol` of atoms removed
        '''
        hydrogens = []
        for atom in self.atoms[:]:
            if atom.symbol != 'H' or len(atom.bonds) != 1:
                continue
            neigh = atom.bond_partners[0]
            if neigh.symbol not in ('C', 'Si'):
                continue
            neigh.mass += atom.mass
            neigh.charge += atom.charge
            hydrogens.append(atom)

        ids_hydrogens = [atom.id_in_mol for atom in hydrogens]
        set_hydrogens = set(hydrogens)

        # cannot modify _atoms and _bonds directly because there are internal data structures to maintain
        # have to call remove_atoms()
        # remove_atoms() will also remove relevant bonds
        self.remove_atoms(hydrogens, update_topology=False)

        # it's faster to reconstruct angle list than iteratively calling _angles.remove()
        self._angles = [conn for conn in self.angles if not set(conn.atoms).intersection(set_hydrogens)]
        self._dihedrals = [conn for conn in self.dihedrals if not set(conn.atoms).intersection(set_hydrogens)]
        self._impropers = [conn for conn in self.impropers if not set(conn.atoms).intersection(set_hydrogens)]

        if self._topology is not None and update_topology:
            self._topology.update_molecules(self._topology.molecules, deepcopy=False)

        return ids_hydrogens

    def add_residue(self, name, atoms, refresh_residues=True):
        '''
        Put a group of atoms into a new residue. These atoms will be removed from their old residues.

        Make sure that these atoms belong to this molecule.
        For performance issue, this is not checked.

        Parameters
        ----------
        name : str
        atoms : list of Atom
        refresh_residues: bool
            If True, the residue list of this molecule will be updated. The `id` and `id_in_mol` for all residues will be assigned.
            This operation can be slow. If you have a lot of residues to add, it is more efficient to set it to False,
            and call `refresh_residues` manually after all residues are added.

        Returns
        -------
        residue : Residue
        '''
        residue = Residue(name)
        for atom in atoms:
            atom.residue._remove_atom(atom)
            residue._add_atom(atom)
        self._added_residues.append(residue)
        if refresh_residues:
            self.refresh_residues()

        return residue

    def remove_residue(self, residue, refresh_residues=True):
        '''
        Remove a residue from this molecule, and put the relevant atoms into the default residue

        Make sure that this residue belongs to this molecule.
        For performance issue, this is not checked.

        Parameters
        ----------
        residue : Residue
        refresh_residues: bool
            If True, the residue list of this molecule will be updated. The `id` and `id_in_mol` for all residues will be assigned.
            This operation can be slow. If you have a lot of residues to remove, it is more efficient to set it to False,
            and call `refresh_residues` manually after all residues are removed.
        '''
        for atom in residue.atoms[:]:
            self._default_residue._add_atom(atom)
        self._added_residues.remove(residue)
        if refresh_residues:
            self.refresh_residues()

    def refresh_residues(self, update_topology=True):
        '''
        Remove empty residues, update `id_in_mol` attributes of each residue in this molecule

        Parameters
        ----------
        update_topology : bool
            If True, the topology this molecule belongs to will assign `id` for all residues
        '''
        self._added_residues = [res for res in self._added_residues if res.n_atom != 0]
        for i, residue in enumerate(self.residues):
            residue.id_in_mol = i
        if self._topology is not None and update_topology:
            for i, residue in enumerate(self._topology.residues):
                residue.id = i

    def add_bond(self, atom1, atom2, order=Bond.Order.UNSPECIFIED, check_existence=False):
        '''
        Add a bond between two atoms.

        Make sure that both these two atoms belong to this molecule.
        For performance issue, this is not checked.

        Parameters
        ----------
        atom1 : Atom
        atom2 : Atom
        check_existence : bool
            If set to True and there is already bond between these two atoms, then do nothing and return None

        Returns
        -------
        bond : [Bond, None]
        '''
        bond = Bond(atom1, atom2, order)
        if check_existence and any(b.equals(bond) for b in self._bonds):
            return None

        self._bonds.append(bond)
        atom1._bonds.append(bond)
        atom2._bonds.append(bond)
        self._is_rdmol_valid = False

        return bond

    def add_angle(self, atom1, atom2, atom3, check_existence=False):
        '''
        Add a angle between three atoms.

        The second atom is the central atom.
        Make sure that both these three atoms belong to this molecule.
        For performance issue, this is not checked.

        Parameters
        ----------
        atom1 : Atom
        atom2 : Atom
        atom3 : Atom
        check_existence : bool
            If set to True and there is already angle between these three atoms, then do nothing and return None
        Returns
        -------
        angle : [Angle, None]
        '''
        angle = Angle(atom1, atom2, atom3)
        if check_existence and any(a.equals(angle) for a in self._angles):
            return None

        self._angles.append(angle)
        return angle

    def add_dihedral(self, atom1, atom2, atom3, atom4, check_existence=False):
        '''
        Add a dihedral between four atoms.

        Make sure that both these four atoms belong to this molecule.
        For performance issue, this is not checked.

        Parameters
        ----------
        atom1 : Atom
        atom2 : Atom
        atom3 : Atom
        atom4 : Atom
        check_existence : bool
            If set to True and there is already dihedral between these three atoms, then do nothing and return None
        Returns
        -------
        dihedral : [Dihedral, None]
        '''
        dihedral = Dihedral(atom1, atom2, atom3, atom4)
        if check_existence and any(d.equals(dihedral) for d in self._dihedrals):
            return None

        self._dihedrals.append(dihedral)
        return dihedral

    def add_improper(self, atom1, atom2, atom3, atom4, check_existence=False):
        '''
        Add a improper between four atoms.

        The fist atom is the central atom.
        Make sure that both these four atoms belong to this molecule.
        For performance issue, this is not checked.

        Parameters
        ----------
        atom1 : Atom
        atom2 : Atom
        atom3 : Atom
        atom4 : Atom
        check_existence : bool
            If set to True and there is already improper between these three atoms, then do nothing and return None
        Returns
        -------
        dihedral : [Improper, None]
        '''
        improper = Improper(atom1, atom2, atom3, atom4)
        if check_existence and any(i.equals(improper) for i in self._impropers):
            return None

        self._impropers.append(improper)
        return improper

    def remove_connectivity(self, connectivity):
        '''
        Remove a connectivity (bond, angle, diheral or improper) from this molecule.

        Make sure that this connectivity belongs to this molecule.
        For performance issue, this is not checked.

        Note that when a bond get removed, the relevant angles, dihedrals and impropers are still there.
        Usually you have to call `generate_angle_dihedral_improper` to keep consistency of connectivities.

        # TODO This operation is slow. A batch version is required for better performance

        Parameters
        ----------
        connectivity : [Bond, Angle, Dihedral, Improper]
        '''
        if type(connectivity) is Bond:
            bond = connectivity
            self._bonds.remove(bond)
            bond.atom1._bonds.remove(bond)
            bond.atom2._bonds.remove(bond)
            self._is_rdmol_valid = False
        elif type(connectivity) is Angle:
            self._angles.remove(connectivity)
        elif type(connectivity) is Dihedral:
            self._dihedrals.remove(connectivity)
        elif type(connectivity) is Improper:
            self._impropers.remove(connectivity)
        else:
            raise Exception('Invalid connectivity')

    def is_similar_to(self, other):
        '''
        Check if this molecule is similar to another molecule.

        It requires two molecules contains the same number of atoms.
        The correspond atoms should have same atom symbol, type and charge.
        The bonds should also be the same.
        But it doesn't consider angles, dihedrals and impropers.

        Parameters
        ----------
        other : Molecule

        Returns
        -------
        is : bool
        '''
        other: Molecule
        if self.n_atom != other.n_atom:
            return False
        if self.n_bond != other.n_bond:
            return False
        for i in range(self.n_atom):
            atom1 = self._atoms[i]
            atom2 = other._atoms[i]
            if atom1.symbol != atom2.symbol or atom1.type != atom2.type or atom1.charge != atom2.charge:
                return False
            if len(atom1._bonds) != len(atom2._bonds):
                return False
            if set(p.id_in_mol for p in atom1.bond_partners) != \
                    set(p.id_in_mol for p in atom2.bond_partners):
                return False
        return True

    def get_adjacency_matrix(self):
        matrix = np.zeros([self.n_atom, self.n_atom], dtype=bool)
        for bond in self.bonds:
            a1, a2 = bond.atom1.id_in_mol, bond.atom2.id_in_mol
            matrix[a1][a2] = True
            matrix[a2][a1] = True

        # virtual site and parent atom are considered as adjacent even though there's no bond between them
        for parent, vsite_atom in self.get_virtual_site_pairs():
            a1, a2 = parent.id_in_mol, vsite_atom.id_in_mol
            matrix[a1][a2] = True
            matrix[a2][a1] = True

        return matrix

    def get_distance_matrix(self, max_bond=None):
        connections = [set() for _ in range(self.n_atom)]
        for bond in self._bonds:
            a1, a2 = bond.atom1.id_in_mol, bond.atom2.id_in_mol
            connections[a1].add(a2)
            connections[a2].add(a1)
        mat = np.zeros((self.n_atom, self.n_atom), dtype=int)

        def fill_matrix(level, center, neighbors, flags):
            if max_bond and level > max_bond:
                return

            if all(flags):
                return

            for j in neighbors:
                if not flags[j]:
                    mat[center][j] = level
                    mat[j][center] = level
                    flags[j] = True

            neighbors_deeper = set()
            for j in neighbors:
                neighbors_deeper.update(connections[j])
            fill_matrix(level + 1, center, neighbors_deeper, flags)

        for i in range(self.n_atom):
            fill_matrix(1, i, connections[i], [j == i for j in range(self.n_atom)])

        return mat

    @property
    def n_atom(self):
        '''
        Number of atoms belong to this molecule

        Returns
        -------
        n : int
        '''
        return len(self._atoms)

    @property
    def n_bond(self):
        '''
        Number of bonds belong to this molecule

        Returns
        -------
        n : int
        '''
        return len(self._bonds)

    @property
    def n_angle(self):
        '''
        Number of angles belong to this molecule

        Returns
        -------
        n : int
        '''
        return len(self._angles)

    @property
    def n_dihedral(self):
        '''
        Number of dihedrals belong to this molecule

        Returns
        -------
        n : int
        '''
        return len(self._dihedrals)

    @property
    def n_improper(self):
        '''
        Number of impropers belong to this molecule

        Returns
        -------
        n : int
        '''
        return len(self._impropers)

    @property
    def n_residue(self):
        return len(self.residues)

    @property
    def atoms(self):
        '''
        List of atoms belong to this molecule

        Returns
        -------
        atoms: list of Atom
        '''
        return self._atoms

    @property
    def bonds(self):
        '''
        List of bonds belong to this molecule

        Returns
        -------
        bonds : list of Bond
        '''
        return self._bonds

    @property
    def angles(self):
        '''
        List of angles belong to this molecule

        Returns
        -------
        angles : list of Angle
        '''
        return self._angles

    @property
    def dihedrals(self):
        '''
        List of dihedrals belong to this molecule

        Returns
        -------
        dihedrals : list of Dihedral
        '''
        return self._dihedrals

    @property
    def impropers(self):
        '''
        List of impropers belong to this molecule

        Returns
        -------
        impropers : list of Improper
        '''
        return self._impropers

    @property
    def residues(self):
        '''
        All the residues in this molecule

        Returns
        -------
        residues : list of Residue
        '''
        if len(self._default_residue.atoms) == 0:
            return self._added_residues
        else:
            return [self._default_residue] + self._added_residues

    @property
    def has_position(self):
        '''
        Whether or not all the atoms in the molecule have positions

        Returns
        -------
        has : bool
        '''
        return all(atom.has_position for atom in self.atoms)

    @property
    def positions(self):
        '''
        Positions of all the atoms in this molecule

        Returns
        -------
        positions : array_like
        '''
        return np.array([atom.position for atom in self.atoms])

    def set_positions(self, positions):
        '''
        Set the positions of all atoms in this molecule

        Parameters
        ----------
        positions : array_like
        '''
        if self.n_atom != len(positions):
            raise Exception('Length of positions should equal to the number of atoms')
        for i, atom in enumerate(self.atoms):
            atom.position = positions[i]

    def get_drude_pairs(self):
        '''
        Retrieve all the Drude dipole pairs belong to this molecule

        Returns
        -------
        pairs :  list of tuple of Atom
            [(parent, drude)]
        '''
        pairs = []
        for bond in self._bonds:
            if bond.atom1.is_drude:
                pairs.append((bond.atom2, bond.atom1))
            elif bond.atom2.is_drude:
                pairs.append((bond.atom1, bond.atom2))
        return pairs

    def get_virtual_site_pairs(self):
        '''
        Retrieve all the virtual site pairs belong to this molecule

        Returns
        -------
        pairs :  list of tuple of Atom
            [(parent, atom_virtual_site)]
        '''
        pairs = []
        for atom in self._atoms:
            if atom.virtual_site is not None:
                pairs.append((atom.virtual_site.parents[0], atom))
        return pairs

    def get_12_13_14_pairs(self):
        '''
        Retrieve all the 1-2, 1-3 and 1-4 pairs based on the bond information.

        The pairs only concerns real atoms. Drude particles  will be ignored.

        Returns
        -------
        pairs12 : list of tuple of Atom
        pairs13 : list of tuple of Atom
        pairs14 : list of tuple of Atom
        '''
        pair_12_set = set()
        pair_13_set = set()
        pair_14_set = set()

        for atom in [a for a in self._atoms if not a.is_drude]:
            partners = [p for p in atom.bond_partners if not p.is_drude]
            for a1, a3 in itertools.combinations(partners, 2):
                pair = tuple(sorted([a1, a3]))
                pair_13_set.add(pair)

        for bond in filter(lambda x: not x.is_drude, self.bonds):
            a2, a3 = bond.atom1, bond.atom2
            pair = tuple(sorted([a2, a3]))
            pair_12_set.add(pair)
            for a1 in [p for p in a2.bond_partners if not p.is_drude]:
                for a4 in [p for p in a3.bond_partners if not p.is_drude]:
                    if a1 != a3 and a2 != a4 and a1 != a4:
                        pair = tuple(sorted([a1, a4]))
                        pair_14_set.add(pair)

        pair_12_list = list(sorted(pair_12_set))
        pair_13_list = list(sorted(pair_13_set - pair_12_set))
        pair_14_list = list(sorted(pair_14_set - pair_13_set.union(pair_12_set)))
        return pair_12_list, pair_13_list, pair_14_list

    def generate_angle_dihedral_improper(self, dihedral=True, improper=True):
        '''
        Generate angle, dihedral and improper from bonds
        The existing angles, dihedrals and impropers will be removed first
        The atoms and bonds concerning Drude particles will be ignored

        Parameters
        ----------
        dihedral: bool
            Whether or not generate dihedrals based on bonds
        improper: bool
            Whether or not generate impropers based on bonds
        '''
        self._angles = []
        self._dihedrals = []
        self._impropers = []

        for atom in [a for a in self._atoms if not a.is_drude]:
            partners = [p for p in atom.bond_partners if not p.is_drude]
            for p1, p2 in itertools.combinations(partners, 2):
                self.add_angle(p1, atom, p2)
            if improper and len(partners) == 3:
                self.add_improper(atom, *sorted(partners))

        if dihedral:
            for bond in filter(lambda x: not x.is_drude, self._bonds):
                atom2 = bond.atom1
                atom3 = bond.atom2
                partners2 = [p for p in atom2.bond_partners if not p.is_drude]
                partners3 = [p for p in atom3.bond_partners if not p.is_drude]
                for atom1, atom4 in itertools.product(partners2, partners3):
                    if atom1 != atom3 and atom2 != atom4 and atom1 != atom4:
                        self.add_dihedral(atom1, atom2, atom3, atom4)

    def guess_connectivity_from_ff(self, ff, bond_limit=0.25, bond_tolerance=0.025, angle_tolerance=None,
                                   pbc='', cell=None):
        '''
        Guess bonds, angles, dihedrals and impropers from force field.

        It requires that atoms types are defined and positions are available.
        The distance between nearby atoms will be calculated.
        If it's smaller than bond_length_limit, then it will be compared with the equilibrium length in FF.
        The bond will be added if a BondTerm is found in FF and the deviation is smaller than bond_tolerance.
        Then angles will be constructed from bonds. If angle_tolerance is None, all angles will be added.
        If angle_tolerance is set (as degree), then AngleTerm must be provided for these angles.
        The angle will be added only if the deviation between angle and equilibrium value in FF is smaller than angle_tolerance.
        Dihedrals and impropers will be constructed form bonds and be added if relevant terms are presented in FF.

        PBC is supported for determining bonds across the periodic cell
        This is useful for simulating infinite structures
        pbc can be '', 'x', 'y', 'xy', 'xz', 'xyz', which means check bonds cross specific boundaries
        cell should also be provided if pbc is not ''

        TODO Add support for triclinic cell

        Parameters
        ----------
        ff : ForceField
        bond_limit : float
        bond_tolerance : float
        angle_tolerance : float
        pbc : str
        cell : UnitCell
        '''
        if not self.has_position:
            raise Exception('Positions are required for guessing connectivity')

        if any(atom.is_drude for atom in self._atoms):
            raise Exception('Drude particles should be removed before guess connectivity')

        if pbc != '':
            if cell is None or cell.volume == 0:
                raise Exception('PBC required but valid cell not provided')
            elif not cell.is_rectangular:
                raise Exception('Triclinic cell haven\'t been implemented')
            else:
                box = cell.get_size()

        self._bonds = []
        self._angles = []
        self._dihedrals = []
        self._impropers = []

        for i in range(self.n_atom):
            atom1 = self.atoms[i]
            try:
                at1 = ff.atom_types[atom1.type].eqt_bond
            except:
                raise Exception(f'AtomType {atom1.type} not found in FF')

            for j in range(i, self.n_atom):
                atom2 = self.atoms[j]
                try:
                    at2 = ff.atom_types[atom2.type].eqt_bond
                except:
                    raise Exception(f'AtomType {atom2.type} not found in FF')

                delta = atom2.position - atom1.position
                if pbc != '':
                    if any((np.abs(delta) > bond_limit) & (np.abs(delta) < box - bond_limit)):
                        continue
                    if 'x' in pbc:
                        delta[0] -= math.ceil(delta[0] / box[0] - 0.5) * box[0]
                    if 'y' in pbc:
                        delta[1] -= math.ceil(delta[1] / box[1] - 0.5) * box[1]
                    if 'z' in pbc:
                        delta[2] -= math.ceil(delta[2] / box[2] - 0.5) * box[2]
                if any(np.abs(delta) > bond_limit):
                    continue

                bterm = BondTerm(at1, at2, 0)
                if bterm.name not in ff.bond_terms.keys():
                    continue

                bterm: BondTerm = ff.bond_terms[bterm.name]
                if any(delta - bterm.length > bond_tolerance):
                    continue

                if abs(np.sqrt(delta.dot(delta)) - bterm.length) <= bond_tolerance:
                    self.add_bond(atom1, atom2)

        # generate angles etc..., and then remove them if requirements are not satisfied
        self.generate_angle_dihedral_improper()

        angles_removed = []
        dihedrals_removed = []
        impropers_removed = []
        if angle_tolerance is not None:
            for angle in self._angles[:]:
                at1 = ff.atom_types[angle.atom1.type].eqt_ang_s
                at2 = ff.atom_types[angle.atom2.type].eqt_ang_c
                at3 = ff.atom_types[angle.atom3.type].eqt_ang_s
                aterm = AngleTerm(at1, at2, at3, 0)
                if aterm.name not in ff.angle_terms.keys():
                    raise Exception(
                        f'{str(angle)} constructed but {str(aterm)} not found in FF')

                aterm: AngleTerm = ff.angle_terms[aterm.name]
                delta21 = angle.atom1.position - angle.atom2.position
                delta23 = angle.atom3.position - angle.atom2.position
                if 'x' in pbc:
                    delta21[0] -= math.ceil(delta21[0] / box[0] - 0.5) * box[0]
                    delta23[0] -= math.ceil(delta23[0] / box[0] - 0.5) * box[0]
                if 'y' in pbc:
                    delta21[1] -= math.ceil(delta21[1] / box[1] - 0.5) * box[1]
                    delta23[1] -= math.ceil(delta23[1] / box[1] - 0.5) * box[1]
                if 'z' in pbc:
                    delta21[2] -= math.ceil(delta21[2] / box[2] - 0.5) * box[2]
                    delta23[2] -= math.ceil(delta23[2] / box[2] - 0.5) * box[2]
                cos = delta21.dot(delta23) / math.sqrt(delta21.dot(delta21) * delta23.dot(delta23))
                theta = np.arccos(np.clip(cos, -1, 1))
                if abs(theta - aterm.theta) > angle_tolerance * DEG2RAD:
                    self.remove_connectivity(angle)
                    angles_removed.append(angle)

        for dihedral in self._dihedrals[:]:
            # consider wildcards in force field
            ats_list = ff.get_eqt_list_for_dihedral(dihedral)
            for ats in ats_list:
                dterm = DihedralTerm(*ats)
                if dterm.name in ff.dihedral_terms.keys():
                    break
            else:
                self.remove_connectivity(dihedral)
                dihedrals_removed.append(dihedral)

        for improper in self._impropers[:]:
            # consider wildcards in force field
            ats_list = ff.get_eqt_list_for_improper(improper)
            for ats in ats_list:
                iterm = ImproperTerm(*ats)
                if iterm.name in ff.improper_terms.keys():
                    break
            else:
                self.remove_connectivity(improper)
                impropers_removed.append(improper)

        if angles_removed != []:
            msg = '%i angles not added because value far from equilibrium: ' \
                  % len(angles_removed) \
                  + ' '.join([i.name for i in angles_removed[:10]])
            if len(angles_removed) > 10:
                msg += ' and more ...'
            logger.warning(msg)

        if dihedrals_removed != []:
            msg = '%i dihedrals not added because parameters not found in FF: ' \
                  % len(dihedrals_removed) \
                  + ' '.join([i.name for i in dihedrals_removed[:10]])
            if len(dihedrals_removed) > 10:
                msg += ' and more ...'
            logger.warning(msg)

        if impropers_removed != []:
            msg = '%i impropers not added because parameters not found in FF: ' \
                  % len(impropers_removed) \
                  + ' '.join([i.name for i in impropers_removed[:10]])
            if len(impropers_removed) > 10:
                msg += ' and more ...'
            logger.warning(msg)

    def generate_drude_particles(self, ff, type_drude='DP_', seed=1, update_topology=True):
        '''
        Generate Drude particles from DrudePolarTerms in force field.

        The atom types should have been defined already.
        Drude particle will not be generated if DrudePolarTerm for its atom type can not be found in the FF.
        Note that The existing Drude particles will be removed before generating.
        The mass defined in the DrudePolarTerm will be transferred from parent atom to the Drude particle.
        The Drude charge will be calculated from the DrudePolarTerm and transferred from parent atom to the Drude particle.
        Bonds between parent-Drude will be generated and added to the topology.
        If AtomType and VdwTerm for generated Drude particles are not found in FF, these terms will be created and added to the FF.

        Parameters
        ----------
        ff : ForceField
        type_drude : str
        seed : int
        update_topology : bool
        '''
        if len(ff.polar_terms) == 0:
            raise Exception('Polar terms not found in force field')

        np.random.seed(seed)

        self.remove_drude_particles(update_topology=False)

        _atype_not_found = set()
        drude_pairs = {}
        for parent in self._atoms:
            atype = ff.atom_types.get(parent.type)
            if atype is None:
                _atype_not_found.add(parent.type)
                continue

            pterm = ff.polar_terms.get(atype.eqt_polar)
            if pterm is None:
                continue
            if type(pterm) is not DrudePolarTerm:
                raise Exception('Polar terms other than DrudePolarTerm haven\'t been implemented')
            drude = Atom()
            drude.is_drude = True
            drude.type = type_drude
            # add Drude particles after all been generated so the name of them are in sequence
            drude.name = 'DP' + str(parent.id_in_mol + 1)
            drude.symbol = 'DP'
            drude.mass = pterm.mass
            parent.mass -= drude.mass
            n_H = len([atom for atom in parent.bond_partners if atom.symbol == 'H'])
            alpha = pterm.alpha + n_H * pterm.merge_alpha_H
            drude.charge = - pterm.get_charge(alpha)
            parent.charge += pterm.get_charge(alpha)
            # update alpha and thole for Drude parent particle
            parent.alpha = alpha
            parent.thole = pterm.thole

            if parent.has_position:
                # make sure Drude and parent atom do not overlap. max deviation 0.005 nm
                drude.position = parent.position + (np.random.random(3) - 0.5) / 100

            drude_pairs[parent] = drude

        if _atype_not_found != set():
            logger.error('%i atom types not found in FF: %s' % (
                len(_atype_not_found), ' '.join(_atype_not_found)))
            raise Exception(f'Generating Drude particles {str(self)} failed')

        for parent, drude in drude_pairs.items():
            self.add_atom(drude, index=self._atoms.index(parent) + 1, update_topology=False)
            self.add_bond(parent, drude)

        if self._topology is not None and update_topology:
            self._topology.update_molecules(self._topology.molecules, deepcopy=False)

        dtype = ff.atom_types.get(type_drude)
        if dtype is None:
            dtype = AtomType(type_drude)
            ff.add_term(dtype)
            logger.warning(f'AtomType for Drude particle not found in FF. '
                           f'{str(dtype)} is added to the FF')
        vdw = LJ126Term(dtype.eqt_vdw, dtype.eqt_vdw, 0.0, 0.0)
        if ff.vdw_terms.get(vdw.name) is None:
            ff.add_term(vdw)
            logger.warning(f'VdwTerm for Drude particle not found in FF. '
                           f'{str(vdw)} with zero interactions is added to the FF')

        for atom in self._atoms:
            if not atom.is_drude and atom.symbol != 'H' and atom not in drude_pairs:
                logger.warning(f'Not all heavy atoms in {str(self)} carry Drude particles')
                break

    def remove_drude_particles(self, update_topology=True):
        '''
        Remove all Drude particles and bonds belong to Drude particles

        The charges and masses carried by Drude particles will be transferred back to parent atoms

        Parameters
        ----------
        update_topology : bool
        '''
        for parent, drude in self.get_drude_pairs():
            parent.mass += drude.mass
            parent.charge += drude.charge
            self.remove_connectivity(drude._bonds[0])
            self.remove_atom(drude, update_topology=False)
        if self._topology is not None and update_topology:
            self._topology.update_molecules(self._topology.molecules, deepcopy=False)

    def generate_virtual_sites(self, ff, update_topology=True):
        '''
        Generate virtual sites from VirtualSiteTerms in force field.

        The atom types should have been defined already.
        Note that The existing virtual sites will be removed before generating.
        The charge won't be assigned by this method.
        Therefore `assign_charge_from_ff` should be called to assign the charges on virtual sites.

        Currently, only TIP4PSiteTerm has been implemented.

        TODO Support other virtual site terms

        Parameters
        ----------
        ff : ForceField
        update_topology : bool
        '''
        if len(ff.virtual_site_terms) == 0:
            raise Exception('Virtual site terms not found in force field')

        self.remove_virtual_sites(update_topology=False)

        for term in ff.virtual_site_terms.values():
            if type(term) is not TIP4PSiteTerm:
                raise Exception('Virtual sites terms other than TIP4PSiteTerm haven\'t been implemented')

        for term in ff.virtual_site_terms.values():
            for angle in self._angles:
                if angle.atom2.type != term.type_O or angle.atom1.type != term.type_H or angle.atom3.type != term.type_H:
                    continue
                atom_vsite = Atom('VS' + str(angle.atom2.id_in_mol + 1))
                atom_vsite.symbol = 'VS'
                atom_vsite.type = term.type
                atom_vsite.virtual_site = TIP4PSite([angle.atom2, angle.atom1, angle.atom3], [term.d])
                atom_vsite.position = atom_vsite.virtual_site.calc_position()
                self.add_atom(atom_vsite, update_topology=False)

        if self._topology is not None and update_topology:
            self._topology.update_molecules(self._topology.molecules, deepcopy=False)

    def remove_virtual_sites(self, update_topology=True):
        '''
        Remove all virtual sites.

        Parameters
        ----------
        update_topology : bool
        '''
        for atom in self._atoms:
            if atom.virtual_site is not None:
                self.remove_atom(atom, update_topology=False)
        if self._topology is not None and update_topology:
            self._topology.update_molecules(self._topology.molecules, deepcopy=False)

    def get_sub_molecule(self, indexes, deepcopy=True):
        '''
        Extract a substructure from this molecule by indexes of atoms.

        The substructure will not contain any bond, angle, dihedral and improper between atoms in substructure and remaining parts.
        Because of this, this method is inefficient. Avoid using this method when integrity check for connectivities is not necessary.

        Residue information will be reconstructed.

        Parameters
        ----------
        indexes : list of int
            The atoms in the substructure will be in the same order as in indexes
        deepcopy : bool
            If set to False, then the atoms and connections in the substructure will be the identical object as the atoms and connections in this molecule.
            The data structure in this molecule will be messed up, and should not be accessed later.

        Returns
        -------
        substructure : Molecule
        '''
        indexes = list(dict.fromkeys(indexes))

        if deepcopy:
            mol = copy.deepcopy(self)
        else:
            mol = self

        # assign temporary atom id. Use it to determine which connectivity to keep
        # cannot use id_in_mol because sub.add_atom will change id_in_mol
        for i, atom in enumerate(mol.atoms):
            atom.id = i

        # copy atoms and reconstruct residues
        res_atoms = {}  # {residue: [atom]}
        for i in indexes:
            atom = mol.atoms[i]
            if atom.residue not in res_atoms:
                res_atoms[atom.residue] = []
            res_atoms[atom.residue].append(atom)

        sub = Molecule()
        for residue, atoms in res_atoms.items():
            for atom in atoms:
                sub.add_atom(atom, update_topology=False)
            if len(res_atoms) > 1:
                sub.add_residue(residue.name, atoms, refresh_residues=False)
            else:
                sub.name = residue.name
        sub.refresh_residues()

        # copy relevant connectivities
        ids_set = set(indexes)
        bonds_kept = {}
        for atom in sub.atoms:
            for bond in reversed(atom.bonds):
                if {a.id for a in bond.atoms} <= ids_set:
                    bonds_kept[bond] = None  # use dict name to store unique bonds in order
                else:
                    atom._bonds.pop()
        sub._bonds.extend(bonds_kept)

        for conn in mol.angles:
            if {a.id for a in conn.atoms} <= ids_set:
                sub._angles.append(conn)
        for conn in mol.dihedrals:
            if {a.id for a in conn.atoms} <= ids_set:
                sub._dihedrals.append(conn)
        for conn in mol.impropers:
            if {a.id for a in conn.atoms} <= ids_set:
                sub._impropers.append(conn)

        return sub

    @staticmethod
    def merge(molecules, deepcopy=True):
        '''
        Merge several molecules into a single molecule.

        Parameters
        ----------
        molecules : list of Molecule
        deepcopy: bool
            If True, the molecules will be deep-copied, and be intact after the mergence.
            Otherwise, the atoms of the merged molecule and of the original molecules will be the same objects.
            Then the original molecules will be unusable.

        Returns
        -------
        merged : Molecule
        '''
        merged = Molecule()
        for mol in molecules:
            if deepcopy:
                m_copy = copy.deepcopy(mol)
            else:
                m_copy = mol
            # should always call `add_atom()` instead of manipulating `_atoms` directly
            for atom in m_copy.atoms:
                merged.add_atom(atom)
            # `atom._bonds` is not empty. Before calling `merged.add_bond()`, I have to empty `atom._bonds`, otherwise `atom._bonds` will contain duplicate bonds.
            # A more efficient way is to directly append the bonds to `merged._bonds`.
            merged._bonds.extend(m_copy._bonds)
            merged._angles.extend(m_copy._angles)
            merged._dihedrals.extend(m_copy._dihedrals)
            merged._impropers.extend(m_copy._impropers)

            # all atoms go to the default residue after `add_atom`, but the old residue still holds the atom list
            for residue in m_copy.residues:
                merged.add_residue(residue.name, residue.atoms, refresh_residues=False)

        merged.refresh_residues()

        return merged

    def split(self, consecutive=False):
        '''
        Split the molecule into smaller pieces based on bond network.

        The atoms in each piece will preserve the original order.
        However, the atoms at the end of original molecule may end up in a piece in the beginning,
        causing the order of all atoms in all the pieces different from original order.
        To avoid this, set consecutive to True.
        In this case, it will make sure all atoms in front pieces will have atom id smaller than atoms in back pieces.

        Residue information will be reconstructed for each piece.

        Parameters
        ----------
        consecutive : bool

        Returns
        -------
        molecules : list of Molecule
        '''
        links = [set() for _ in range(self.n_atom)]
        for bond in self.bonds:
            id1, id2 = bond.atom1.id_in_mol, bond.atom2.id_in_mol
            links[id1].add(id2)
            links[id2].add(id1)
        if consecutive:
            for i in range(self.n_atom):
                min_neigh = min(links[i] | {i})
                max_neigh = max(links[i] | {i})
                links[i] = set(range(min_neigh, max_neigh + 1))

        clusters = find_clusters_in_graph(self.n_atom, links)
        for cluster in clusters:
            cluster.sort()
        clusters.sort()

        mol = copy.deepcopy(self)

        # A convenient way is to call get_sub_molecule() for each cluster.
        # But get_sub_molecule() is inefficient because it checks the integrity of connectivities.
        # This check is unnecessary since the herein clusters are split based on bonds.
        #
        # pieces = [mol.get_sub_molecule(cluster, deepcopy=False) for cluster in clusters]
        #
        # return pieces

        pieces = []
        for cluster in clusters:
            # copy atoms and reconstruct residues
            res_atoms = {}  # {residue: [atom]}
            for i in cluster:
                atom = mol.atoms[i]
                if atom.residue not in res_atoms:
                    res_atoms[atom.residue] = []
                res_atoms[atom.residue].append(atom)

            piece = Molecule()
            for residue, atoms in res_atoms.items():
                for atom in atoms:
                    piece.add_atom(atom, update_topology=False)
                if len(res_atoms) > 1:
                    piece.add_residue(residue.name, atoms, refresh_residues=False)
                else:
                    piece.name = residue.name
            piece.refresh_residues()

            pieces.append(piece)

        # The atom._bonds already contains correct bonds, but the piece._bonds is empty.
        # Cannot call piece.add_bond() because it will add the bond to atom._bonds again.
        # So I have to add bonds to molecule._bonds manually.
        # This is actually more efficient than piece.add_bond().
        for conn in mol.bonds:
            conn.atom1.molecule._bonds.append(conn)
        for conn in mol.angles:
            conn.atom1.molecule._angles.append(conn)
        for conn in mol.dihedrals:
            conn.atom1.molecule._dihedrals.append(conn)
        for conn in mol.impropers:
            conn.atom1.molecule._impropers.append(conn)

        return pieces

    def split_residues(self, consecutive=False):
        '''
        Split the molecule into smaller pieces based on bond network between residues.

        The atoms belonging to the same residue will always be in the same piece.

        The residues in each piece will preserve the original order.
        However, the residues at the end of original molecule may end up in a piece in the beginning,
        causing the order of all residues in all the pieces different from original order.
        To avoid this, set consecutive to True.
        In this case, it will make sure all residues in front pieces will have atom id smaller than residues in back pieces.

        Residue information will be reconstructed for each piece.

        Parameters
        ----------
        consecutive : bool

        Returns
        -------
        molecules : list of Molecule
        '''
        links = [set() for _ in range(self.n_residue)]
        for bond in self.bonds:
            resid1, resid2 = bond.atom1.residue.id_in_mol, bond.atom2.residue.id_in_mol
            if resid1 != resid2:
                links[resid1].add(resid2)
                links[resid2].add(resid1)
        if consecutive:
            for i in range(self.n_residue):
                min_neigh = min(links[i] | {i})
                max_neigh = max(links[i] | {i})
                links[i] = set(range(min_neigh, max_neigh + 1))

        clusters = find_clusters_in_graph(self.n_residue, links)
        for cluster in clusters:
            cluster.sort()
        clusters.sort()

        pieces = []
        mol = copy.deepcopy(self)
        for cluster in clusters:
            piece = Molecule()
            for rid in cluster:
                residue = mol.residues[rid]
                for atom in residue.atoms:
                    piece.add_atom(atom)  # atom._bonds is not empty
                if len(cluster) > 1:
                    piece.add_residue(residue.name, residue.atoms, refresh_residues=False)
                else:
                    piece.name = residue.name
            piece.refresh_residues()
            pieces.append(piece)

        # The atom._bonds already contains correct bonds, but the piece._bonds is empty.
        # Cannot call piece.add_bond() because it will add the bond to atom._bonds again.
        # So I have to add bonds to molecule._bonds manually.
        # This is actually more efficient than piece.add_bond().
        for conn in mol.bonds:
            conn.atom1.molecule._bonds.append(conn)
        for conn in mol.angles:
            conn.atom1.molecule._angles.append(conn)
        for conn in mol.dihedrals:
            conn.atom1.molecule._dihedrals.append(conn)
        for conn in mol.impropers:
            conn.atom1.molecule._impropers.append(conn)

        return pieces
