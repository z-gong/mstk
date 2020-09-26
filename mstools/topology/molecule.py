import math
import itertools
import copy
import numpy as np
from .atom import Atom
from .virtualsite import *
from .connectivity import *
from .unitcell import UnitCell
from ..forcefield import *
from ..utils import create_mol_from_smiles
from .. import logger


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
        self.name = name
        self._topology = None
        self._atoms: [Atom] = []
        self._bonds: [Bond] = []
        self._angles: [Angle] = []
        self._dihedrals: [Dihedral] = []
        self._impropers: [Improper] = []
        self._obmol = None  # this is for typing based on SMARTS

    def __repr__(self):
        return f'<Molecule: {self.name} {self.id}>'

    def __deepcopy__(self, memodict={}):
        '''
        If there are virtual sites, they will be constructed with new atoms also
        '''
        mol = Molecule()
        mol.id = self.id
        mol.name = self.name
        if self._obmol is not None:
            import openbabel
            mol._obmol = openbabel.OBMol(self._obmol)
        for atom in self._atoms:
            atom_new = copy.deepcopy(atom)
            mol.add_atom(atom_new)
        for bond in self._bonds:
            idx1 = bond.atom1.id_in_mol
            idx2 = bond.atom2.id_in_mol
            mol.add_bond(mol._atoms[idx1], mol._atoms[idx2])
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

        for i, atom in enumerate(self._atoms):
            vsite = atom.virtual_site
            if vsite is not None:
                new_atom = mol.atoms[i]
                new_parents = [mol.atoms[p.id_in_mol] for p in vsite.parents]
                new_atom.virtual_site = VirtualSiteFactory.create(vsite.site_type, new_parents, vsite.parameters)

        return mol

    @staticmethod
    def from_smiles(smiles):
        '''
        Initialize a molecule from SMILES string.

        OpenBabel is used for parsing SMILES. The Hydrogen atoms will be created.
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
        smiles = words[0]
        if len(words) > 1:
            name = words[1]
        else:
            name = None

        py_mol = create_mol_from_smiles(smiles)
        mol = Molecule.from_pybel(py_mol, name)

        return mol

    @staticmethod
    def from_pybel(py_mol, name=None):
        '''
        Initialize a molecule from a `pybel Molecule` object.

        Parameters
        ----------
        py_mol : pybel.Molecule
        name : str
            The name of the molecule. If not provided, the formula will be used as the name.

        Returns
        -------
        molecule : Molecule
        '''
        try:
            import openbabel
        except ImportError:
            raise ImportError('OpenBabel is required for parsing SMILES')

        mol = Molecule()
        for i, a in enumerate(py_mol.atoms):
            atom = Atom()
            element = Element(a.atomicnum)
            atom.name = element.symbol + str(i + 1)
            atom.symbol = element.symbol
            atom.mass = element.mass
            atom.formal_charge = a.formalcharge
            atom.position = [i / 10 for i in a.coords]  # convert A to nm
            mol.add_atom(atom)
        for b in openbabel.OBMolBondIter(py_mol.OBMol):
            atom1 = mol.atoms[b.GetBeginAtomIdx() - 1]
            atom2 = mol.atoms[b.GetEndAtomIdx() - 1]
            mol.add_bond(atom1, atom2)
        mol.generate_angle_dihedral_improper()
        mol._obmol = py_mol.OBMol

        if name is not None:
            mol.name = name
        else:
            mol.name = py_mol.formula

        return mol

    @property
    def obmol(self):
        '''
        The `openbabel.OBMol` object associated with this molecule.

        It is required by ZftTyper typing engine, which performs SMARTS matching on the molecule.
        The obmol attribute will be assigned if the molecule is initialized from SMILES or Pybel Molecule.
        If this information is not available, an Exception will be raised.

        Returns
        -------
        obmol : openbabel.OBMol
        '''
        if self._obmol is not None:
            return self._obmol
        else:
            raise Exception('obmol not assigned for this molecule')

    @property
    def topology(self):
        '''
        The topology this molecule belongs to

        Returns
        -------
        topology : Topology

        '''
        return self._topology

    def add_atom(self, atom, index=None, update_topology=True):
        '''
        Add an atom to this molecule.

        The id_in_mol attribute of all atoms will be updated after insertion.

        Parameters
        ----------
        atom : Atom
        index : int or None
            If index is None, the new atom will be the last atom. Otherwise, it will be inserted in front of index-th atom.
        update_topology : bool
            If update_topology is True, the topology this molecule belongs to will update its atom list and assign id for all atoms.
            Otherwise, you have to re-init the topology manually so that the topological information is correct.
        '''
        atom._molecule = self
        if index is None:
            self._atoms.append(atom)
            atom.id_in_mol = len(self._atoms) - 1
        else:
            self._atoms.insert(index, atom)
            for i, at in enumerate(self._atoms):
                at.id_in_mol = i
        if self._topology is not None and update_topology:
            self._topology.update_molecules(self._topology.molecules, deepcopy=False)

    def remove_atom(self, atom, update_topology=True):
        '''
        Remove an atom and all the bonds connected to the atom from this molecule.
        The id_in_mol attribute of all atoms will be updated after removal.

        Parameters
        ----------
        atom : Atom
        update_topology : bool
            If update_topology is True, the topology this molecule belongs to will update its atom list and assign id for all atoms.
            Otherwise, you have to re-init the topology manually so that the topological information is correct.
        '''
        for bond in atom._bonds[:]:
            self.remove_connectivity(bond)
        self._atoms.remove(atom)
        atom._molecule = None
        for i, at in enumerate(self._atoms):
            at.id_in_mol = i
        if self._topology is not None and update_topology:
            self._topology.update_molecules(self._topology.molecules, deepcopy=False)

    def add_bond(self, atom1, atom2, check_existence=False):
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
        bond = Bond(atom1, atom2)
        if check_existence and bond in self._bonds:
            return None

        self._bonds.append(bond)
        atom1._bonds.append(bond)
        atom2._bonds.append(bond)
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
        if check_existence and angle in self._angles:
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
        if check_existence and dihedral in self._dihedrals:
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
        if check_existence and improper in self._impropers:
            return None

        self._impropers.append(improper)
        return improper

    def remove_connectivity(self, connectivity):
        '''
        Remove a connectivity (bond, angle, diheral or improper) from this molecule.

        Make sure that this connectivity belongs to this molecule.
        For performance issue, this is not checked.

        Parameters
        ----------
        connectivity : [Bond, Angle, Dihedral, Improper]
        '''
        if type(connectivity) is Bond:
            bond = connectivity
            self._bonds.remove(bond)
            bond.atom1._bonds.remove(bond)
            bond.atom2._bonds.remove(bond)
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
    def has_position(self):
        '''
        Whether or not all the atoms in the molecule have positions

        Returns
        -------
        has : bool
        '''
        return all(atom.has_position for atom in self.atoms)

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

    def generate_angle_dihedral_improper(self):
        '''
        Generate angle, dihedral and improper from bonds
        The existing angles, dihedrals and impropers will be removed first
        The atoms and bonds concerning Drude particles will be ignored
        '''
        self._angles = []
        self._dihedrals = []
        self._impropers = []

        for atom in [a for a in self._atoms if not a.is_drude]:
            partners = [p for p in atom.bond_partners if not p.is_drude]
            for p1, p2 in itertools.combinations(partners, 2):
                self.add_angle(p1, atom, p2)
            if len(partners) == 3:
                self.add_improper(atom, *sorted(partners))

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

        * TODO Add support for triclinic cell

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
                box = cell.size

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
                    if any((delta > bond_limit) & (delta < box - bond_limit)):
                        continue
                    if 'x' in pbc:
                        delta[0] -= math.ceil(delta[0] / box[0] - 0.5) * box[0]
                    if 'y' in pbc:
                        delta[1] -= math.ceil(delta[1] / box[1] - 0.5) * box[1]
                    if 'z' in pbc:
                        delta[2] -= math.ceil(delta[2] / box[2] - 0.5) * box[2]
                if any(delta > bond_limit):
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
                if abs(theta * 180 / math.pi - aterm.theta) > angle_tolerance:
                    self.remove_connectivity(angle)
                    angles_removed.append(angle)

        for dihedral in self._dihedrals[:]:
            at1 = ff.atom_types[dihedral.atom1.type].eqt_dih_s
            at2 = ff.atom_types[dihedral.atom2.type].eqt_dih_c
            at3 = ff.atom_types[dihedral.atom3.type].eqt_dih_c
            at4 = ff.atom_types[dihedral.atom4.type].eqt_dih_s
            dterm = DihedralTerm(at1, at2, at3, at4)
            if dterm.name not in ff.dihedral_terms.keys():
                self.remove_connectivity(dihedral)
                dihedrals_removed.append(dihedral)

        for improper in self._impropers[:]:
            at1 = ff.atom_types[improper.atom1.type].eqt_imp_c
            at2 = ff.atom_types[improper.atom2.type].eqt_imp_s
            at3 = ff.atom_types[improper.atom3.type].eqt_imp_s
            at4 = ff.atom_types[improper.atom4.type].eqt_imp_s
            iterm = ImproperTerm(at1, at2, at3, at4)
            if iterm.name not in ff.improper_terms.keys():
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
        Generate Drude particles from DrudeTerms in force field.

        The atom types should have been defined already.
        Drude particle will not be generated if DrudeTerm for its atom type can not be found in the FF.
        Note that The existing Drude particles will be removed before generating.
        The mass defined in the DrudeTerm will be transferred from parent atom to the Drude particle.
        The Drude charge will be calculated from the DrudeTerm and transferred from parent atom to the Drude particle.
        Bonds between parent-Drude will be generated and added to the topology.
        If AtomType and VdwTerm for generated Drude particles are not found in FF, these terms will be created and added to the FF.

        Parameters
        ----------
        ff : ForceField
        type_drude : str
        seed : int
        update_topology : bool
        '''
        if len(ff.polarizable_terms) == 0:
            raise Exception('Polarizable terms not found in force field')

        np.random.seed(seed)

        self.remove_drude_particles(update_topology=False)

        _atype_not_found = set()
        drude_pairs = {}
        for parent in self._atoms:
            atype = ff.atom_types.get(parent.type)
            if atype is None:
                _atype_not_found.add(parent.type)
                continue

            pterm = ff.polarizable_terms.get(atype.eqt_polar)
            if pterm is None:
                continue
            if type(pterm) is not DrudeTerm:
                raise Exception('Polarizable terms other than DrudeTerm haven\'t been implemented')
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

    def assign_mass_from_ff(self, ff):
        '''
        Assign masses of all atoms from the force field.

        The atom types should have been defined, and the AtomType in FF should carry mass information.
        The masses of Drude particles will be determined ONLY from the DrudeTerm,
        and the the same amount of mass will be subtracted from the parent atoms.
        That is, if there is an AtomType for Drude particles, the mass attribute of this AtomType will be simply ignored.

        Parameters
        ----------
        ff : ForceField
        '''
        _atype_not_found = set()
        _atype_no_mass = set()
        _zero_mass = set()
        for atom in self._atoms:
            if atom.is_drude:
                continue
            atype = ff.atom_types.get(atom.type)
            if atype is None:
                _atype_not_found.add(atom.type)
                continue
            elif atype.mass == -1:
                _atype_no_mass.add(atom.type)
            elif atype.mass == 0:
                _zero_mass.add(atom.type)
            atom.mass = atype.mass

        if _atype_not_found != set():
            logger.error('%i atom types not found: %s' % (
                len(_atype_not_found), ' '.join(_atype_not_found)))
            raise Exception(f'Assign mass for {str(self)} failed')
        if _atype_no_mass != set():
            logger.error('%i atom types in FF do not carry mass information: %s' % (
                len(_atype_no_mass), ' '.join(_atype_no_mass)))
            raise Exception(f'Assign mass {str(self)}  failed')
        if _zero_mass != set():
            logger.warning(f'%i atoms in {str(self)} carry zero mass: %s' % (
                len(_zero_mass), ' '.join(_zero_mass)))

        for parent, drude in self.get_drude_pairs():
            atype = ff.atom_types[parent.type]
            pterm = ff.polarizable_terms.get(atype.eqt_polar)
            if pterm is None:
                raise Exception(f'Polarizable term for {atype} not found in FF')
            if type(pterm) != DrudeTerm:
                raise Exception('Polarizable terms other than DrudeTerm haven\'t been implemented')
            drude.mass = pterm.mass
            parent.mass -= pterm.mass

    def assign_charge_from_ff(self, ff, transfer_bci_terms=False):
        '''
        Assign charges of all atoms from the force field.

        The atom types should have been defined.
        The charge of corresponding AtomType in FF will be assigned to the atoms first.
        Then if there are ChargeIncrementTerm, the charge will be transferred between bonded atoms.
        The charge of Drude particles will be determined ONLY from the DrudePolarizableTerm,
        and the the same amount of charge will be subtracted from the parent atoms.
        That is, if there is an AtomType for Drude particles, the charge attribute of this AtomType will be simply ignored.

        If charge increment parameters required by the topology are not found in the force field,
        `mstools` can try to transfer parameters from more general terms in the force field.
        e.g. if bci term between 'c_4o2' and 'n_3' is not found, the bci term between 'c_4' and 'n_3' will be used.
        This mechanism is initially designed for TEAM force field.
        It is disabled by default, and can be enabled by setting argument `transfer_bci_terms` to True.

        Parameters
        ----------
        ff : ForceField
        transfer_bci_terms : bool, optional
        '''
        _q_zero = []
        for atom in self._atoms:
            if atom.is_drude:
                continue
            try:
                charge = ff.atom_types[atom.type].charge
            except:
                raise Exception(f'Atom type {atom.type} not found in FF')
            atom.charge = charge
            if charge == 0:
                _q_zero.append(atom)

        if len(ff.bci_terms) == 0 and len(_q_zero) > 0:
            msg = '%i atoms in %s have zero charge. But charge increment terms not found in FF: %s' % (
                len(_q_zero), self, ' '.join([a.name for a in _q_zero[:10]]))
            if len(_q_zero) > 10:
                msg += ' and more ...'
            logger.warning(msg)

        if len(ff.bci_terms) > 0:
            _qterm_not_found = set()
            _qterm_transferred = set()
            for bond in filter(lambda x: not x.is_drude, self._bonds):
                _found = False
                ats = ff.get_eqt_for_bci(bond)
                if ats[0] == ats[1]:
                    continue
                _term = ChargeIncrementTerm(*ats, 0)
                try:
                    qterm, direction = ff.get_bci_term(*ats)
                    _found = True
                except FFTermNotFoundError:
                    if transfer_bci_terms:
                        qterm, score = dff_fuzzy_match(_term, ff)
                        if qterm is not None:
                            direction = 1 if _term.type1 == ats[0] else -1
                            _found = True
                            _qterm_transferred.add((_term.name, qterm.name, score))

                if _found:
                    increment = qterm.value * direction
                    bond.atom1.charge += increment
                    bond.atom2.charge -= increment
                else:
                    _qterm_not_found.add(_term.name)

            if len(_qterm_transferred) > 0:
                logger.warning('%i charge increment terms not found in FF. Transferred from similar terms with score\n'
                               '        %s' % (len(_qterm_transferred),
                                               '\n        '.join(['%-16s %-16s %i' % x for x in _qterm_transferred])))
            if len(_qterm_not_found) > 0:
                logger.warning('%i charge increment terms not found in FF. Make sure FF is correct\n'
                               '        %s' % (len(_qterm_not_found),
                                               '\n        '.join(_qterm_not_found)))

        for parent, drude in self.get_drude_pairs():
            atype = ff.atom_types[parent.type]
            pterm = ff.polarizable_terms.get(atype.eqt_polar)
            if pterm is None:
                raise Exception(f'Polarizable term for {atype} not found in FF')
            if type(pterm) != DrudeTerm:
                raise Exception('Polarizable terms other than DrudeTerm haven\'t been implemented')
            n_H = len([atom for atom in parent.bond_partners if atom.symbol == 'H'])
            alpha = pterm.alpha + n_H * pterm.merge_alpha_H
            drude.charge = -pterm.get_charge(alpha)
            parent.charge += pterm.get_charge(alpha)
            # update alpha and thole for Drude parent particle
            parent.alpha = alpha
            parent.thole = pterm.thole
