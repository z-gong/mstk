import math
import itertools
import random
import copy
import numpy as np
from .atom import Atom
from .connectivity import *
from .unitcell import UnitCell
from ..forcefield import ForceField, Element
from ..forcefield.ffterm import *
from ..forcefield.errors import *
from .. import logger


class Molecule():
    def __init__(self, name: str = 'UNK'):
        self.id: int = -1  # id in topology. -1 means information haven\'t been initiated
        self.name: str = name
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
            idx1 = bond.atom1._id_in_molecule
            idx2 = bond.atom2._id_in_molecule
            mol.add_bond(mol._atoms[idx1], mol._atoms[idx2])
        for angle in self._angles:
            idx1 = angle.atom1._id_in_molecule
            idx2 = angle.atom2._id_in_molecule
            idx3 = angle.atom3._id_in_molecule
            mol.add_angle(mol._atoms[idx1], mol._atoms[idx2], mol._atoms[idx3])
        for dihedral in self._dihedrals:
            idx1 = dihedral.atom1._id_in_molecule
            idx2 = dihedral.atom2._id_in_molecule
            idx3 = dihedral.atom3._id_in_molecule
            idx4 = dihedral.atom4._id_in_molecule
            mol.add_dihedral(mol._atoms[idx1], mol._atoms[idx2], mol._atoms[idx3], mol._atoms[idx4])
        for improper in self._impropers:
            idx1 = improper.atom1._id_in_molecule
            idx2 = improper.atom2._id_in_molecule
            idx3 = improper.atom3._id_in_molecule
            idx4 = improper.atom4._id_in_molecule
            mol.add_improper(mol._atoms[idx1], mol._atoms[idx2], mol._atoms[idx3], mol._atoms[idx4])
        return mol

    @staticmethod
    def from_smiles(smiles):
        try:
            import pybel
        except ImportError:
            raise ImportError('OpenBabel is required for parsing SMILES')

        words = smiles.strip().split()
        smiles = words[0]
        if len(words) > 1:
            name = words[1]
        else:
            name = None

        try:
            py_mol = pybel.readstring('smi', smiles)
        except:
            raise Exception('Invalid SMILES')

        py_mol.addh()
        py_mol.make3D()
        mol = Molecule.from_pybel(py_mol, name)

        return mol

    @staticmethod
    def from_pybel(py_mol, name=None):
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
            atom.has_position = True
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
        if self._obmol is not None:
            return self._obmol
        else:
            raise Exception('obmol not assigned for this molecule')

    @property
    def topology(self):
        return self._topology

    def add_atom(self, atom: Atom, index=None, update_topology=True):
        '''
        Add an atom to this molecule.
        If index is None, the new atom will be the last atom. Otherwise, it will be inserted in front of index-th atom.
        The _id_in_molecule attribute of all atoms will be updated after insertion.
        If update_topology is True, the topology this molecule belongs to will update its atom list and assign id for all atoms.
        Otherwise, you have to re-init the topology manually so that the topological information is correct.
        '''
        atom._molecule = self
        if index is None:
            self._atoms.append(atom)
            atom._id_in_molecule = len(self._atoms) - 1
        else:
            self._atoms.insert(index, atom)
            for i, at in enumerate(self._atoms):
                at._id_in_molecule = i
        if self._topology is not None and update_topology:
            self._topology.update_molecules(self._topology.molecules, deepcopy=False)

    def remove_atom(self, atom: Atom, update_topology=True):
        '''
        Remove an atom and all the bonds connected to the atom from this molecule.
        The _id_in_molecule attribute of all atoms will be updated after removal.
        If update_topology is True, the topology this molecule belongs to will update its atom list and assign id for all atoms.
        Otherwise, you have to re-init the topology manually so that the topological information is correct.
        '''
        for bond in atom._bonds[:]:
            self.remove_bond(bond)
        self._atoms.remove(atom)
        atom._molecule = None
        for i, at in enumerate(self._atoms):
            at._id_in_molecule = i
        if self._topology is not None and update_topology:
            self._topology.update_molecules(self._topology.molecules, deepcopy=False)

    def add_bond(self, atom1: Atom, atom2: Atom, check_existence=False):
        bond = Bond(atom1, atom2)
        if check_existence and bond in self._bonds:
            return

        self._bonds.append(bond)
        atom1._bonds.append(bond)
        atom2._bonds.append(bond)
        return bond

    def remove_bond(self, bond: Bond):
        self._bonds.remove(bond)
        bond.atom1._bonds.remove(bond)
        bond.atom2._bonds.remove(bond)

    def add_angle(self, atom1: Atom, atom2: Atom, atom3: Atom, check_existence=False):
        angle = Angle(atom1, atom2, atom3)
        if check_existence and angle in self._angles:
            return

        self._angles.append(angle)
        return angle

    def remove_angle(self, angle: Angle):
        self._angles.remove(angle)

    def add_dihedral(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom,
                     check_existence=False):
        dihedral = Dihedral(atom1, atom2, atom3, atom4)
        if check_existence and dihedral in self._dihedrals:
            return

        self._dihedrals.append(dihedral)
        return dihedral

    def remove_dihedral(self, dihedral: Dihedral):
        self._dihedrals.remove(dihedral)

    def add_improper(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom,
                     check_existence=False):
        improper = Improper(atom1, atom2, atom3, atom4)
        if check_existence and improper in self._impropers:
            return

        self._impropers.append(improper)
        return improper

    def remove_improper(self, improper: Improper):
        self._impropers.remove(improper)

    def similar_to(self, other) -> bool:
        '''
        See if this molecule is similar to other molecule.
        It requires two molecules contains the same number of atoms.
        The correspond atoms should have same atom symbol, type and charge.
        The bonds should also be the same.
        But it doesn't consider angles, dihedrals and impropers.
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
            if set(p._id_in_molecule for p in atom1.bond_partners) != \
                    set(p._id_in_molecule for p in atom2.bond_partners):
                return False
        return True

    @property
    def n_atom(self):
        return len(self._atoms)

    @property
    def n_bond(self):
        return len(self._bonds)

    @property
    def n_angle(self):
        return len(self._angles)

    @property
    def n_dihedral(self):
        return len(self._dihedrals)

    @property
    def n_improper(self):
        return len(self._impropers)

    @property
    def atoms(self):
        return self._atoms

    @property
    def bonds(self):
        return self._bonds

    @property
    def angles(self):
        return self._angles

    @property
    def dihedrals(self):
        return self._dihedrals

    @property
    def impropers(self):
        return self._impropers

    @property
    def has_position(self):
        return all(atom.has_position for atom in self.atoms)

    def get_drude_pairs(self) -> [(Atom, Atom)]:
        '''
        [(parent, drude)]
        '''
        pairs = []
        for bond in self._bonds:
            if bond.atom1.is_drude:
                pairs.append((bond.atom2, bond.atom1))
            elif bond.atom2.is_drude:
                pairs.append((bond.atom1, bond.atom2))
        return pairs

    def get_12_13_14_pairs(self) -> ([(Atom, Atom)], [(Atom, Atom)], [(Atom, Atom)]):
        '''
        The pairs only concerns real atoms. Drude particles  will be ignored.
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

    def guess_connectivity_from_ff(self, ff: ForceField, bond_limit=0.25,
                                   bond_tolerance=0.025, angle_tolerance=None,
                                   pbc='', cell: UnitCell = None):
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
                    self.remove_angle(angle)
                    angles_removed.append(angle)

        for dihedral in self._dihedrals[:]:
            at1 = ff.atom_types[dihedral.atom1.type].eqt_dih_s
            at2 = ff.atom_types[dihedral.atom2.type].eqt_dih_c
            at3 = ff.atom_types[dihedral.atom3.type].eqt_dih_c
            at4 = ff.atom_types[dihedral.atom4.type].eqt_dih_s
            dterm = DihedralTerm(at1, at2, at3, at4)
            if dterm.name not in ff.dihedral_terms.keys():
                self.remove_dihedral(dihedral)
                dihedrals_removed.append(dihedral)

        for improper in self._impropers[:]:
            at1 = ff.atom_types[improper.atom1.type].eqt_imp_c
            at2 = ff.atom_types[improper.atom2.type].eqt_imp_s
            at3 = ff.atom_types[improper.atom3.type].eqt_imp_s
            at4 = ff.atom_types[improper.atom4.type].eqt_imp_s
            iterm = ImproperTerm(at1, at2, at3, at4)
            if iterm.name not in ff.improper_terms.keys():
                self.remove_improper(improper)
                impropers_removed.append(improper)

        if angles_removed != []:
            string = '%i angles not added because value far from equilibrium: ' \
                     % len(angles_removed) \
                     + ' '.join([str(i) for i in angles_removed[:10]])
            if len(angles_removed) > 10:
                string += ' and more ...'
            logger.warning(string)

        if dihedrals_removed != []:
            string = '%i dihedrals not added because parameters not found in FF: ' \
                     % len(dihedrals_removed) \
                     + ' '.join([str(i) for i in dihedrals_removed[:10]])
            if len(dihedrals_removed) > 10:
                string += ' and more ...'
            logger.warning(string)

        if impropers_removed != []:
            string = '%i impropers not added because parameters not found in FF: ' \
                     % len(impropers_removed) \
                     + ' '.join([str(i) for i in impropers_removed[:10]])
            if len(impropers_removed) > 10:
                string += ' and more ...'
            logger.warning(string)

    def generate_drude_particles(self, ff: ForceField, type_drude='DP_', seed=1):
        '''
        Generate Drude particles from DrudeTerms in force field.
        The atom types should have been defined already.
        Drude particle will not be generated if DrudeTerm for its atom type can not be found in the FF.
        Note that The existing Drude particles will be removed before generating.
        The mass defined in the DrudeTerm will be transferred from parent atom to the Drude particle.
        The Drude charge will be calculated from the DrudeTerm and transferred from parent atom to the Drude particle.
        Bonds between parent-Drude will be generated and added to the topology.
        If AtomType and VdwTerm for generated Drude particles are not found in FF, these terms will be created and added to the FF.
        '''
        if len(ff.polarizable_terms) == 0:
            raise Exception('Polarizable terms not found in force field')

        random.seed(seed)

        self.remove_drude_particles()
        drude_pairs = []
        for parent in self._atoms[:]:
            atype = ff.atom_types.get(parent.type)
            if atype is None:
                raise Exception(f'Atom type {parent.type} not found in FF')

            pterm = ff.polarizable_terms.get(atype.eqt_polar)
            if pterm is None:
                continue
            if type(pterm) != DrudeTerm:
                raise Exception('Polarizable terms other than DrudeTerm haven\'t been implemented')
            drude = Atom()
            drude.is_drude = True
            drude.type = type_drude
            # add Drude particles after all been generated so the name of them are in sequence
            drude.name = 'DP' + str(parent._id_in_molecule + 1)
            drude.symbol = 'DP'
            drude.mass = pterm.mass
            parent.mass -= drude.mass
            n_H = len([atom for atom in parent.bond_partners if atom.symbol == 'H'])
            alpha = pterm.alpha + n_H * pterm.merge_alpha_H
            drude.charge = - pterm.get_charge(alpha)
            parent.charge += pterm.get_charge(alpha)
            # store _alpha and _thole for PSF exporting
            parent._alpha = alpha
            parent._thole = pterm.thole

            if parent.has_position:
                # make sure Drude and parent atom do not overlap. max deviation 0.001 nm
                deviation = random.random() / 1000
                drude.position = parent.position + [deviation, deviation, deviation]
                drude.has_position = True

            drude_pairs.append((parent, drude))

        for parent, drude in drude_pairs:
            self.add_atom(drude, index=self._atoms.index(parent) + 1, update_topology=False)
            self.add_bond(parent, drude)

        if self._topology is not None:
            self._topology.update_molecules(self._topology.molecules, deepcopy=False)

        dtype = ff.atom_types.get(type_drude)
        if dtype is None:
            dtype = AtomType(type_drude)
            ff.add_term(dtype)
            logger.warning(f'AtomType for Drude particle not found in FF. '
                           f'{str(dtype)} is added to the FF')
        vdw = LJ126Term(dtype.eqt_vdw, dtype.eqt_vdw, 0.0, 0.0)
        vdw = ff.vdw_terms.get(vdw.name)
        if vdw is None:
            ff.add_term(vdw)
            logger.warning(f'VdwTerm for Drude particle not found in FF. '
                           f'{str(vdw)} with zero interactions is added to the FF')

    def remove_drude_particles(self):
        '''
        Remove all Drude particles and bonds belong to Drude particles
        The charges and masses carried by Drude particles will be transferred back to parent atoms
        :return:
        '''
        for parent, drude in self.get_drude_pairs():
            parent.mass += drude.mass
            parent.charge += drude.charge
            self.remove_bond(drude._bonds[0])
            self.remove_atom(drude, update_topology=False)
        if self._topology is not None:
            self._topology.update_molecules(self._topology.molecules, deepcopy=False)

    def assign_mass_from_ff(self, ff: ForceField):
        '''
        Assign masses of all atoms from the force field.
        The atom types should have been defined, and the AtomType in FF should carry mass information.
        The masses of Drude particles will be determined ONLY from the DrudeTerm,
        and the the same amount of mass will be subtracted from the parent atoms.
        That is, if there is an AtomType for Drude particles, the mass attribute of this AtomType will be simply ignored.
        '''
        for atom in self._atoms:
            if atom.is_drude:
                continue
            atype = ff.atom_types.get(atom.type)
            if atype is None:
                raise Exception(f'Atom type {atom.type} not found in FF')
            if atype.mass == -1:
                raise Exception(f'{atype} does not carry mass information')
            if atype.mass == 0:
                logger.warning(f'{atype} carries zero mass. You should make sure it is correct')
            atom.mass = atype.mass

        for parent, drude in self.get_drude_pairs():
            atype = ff.atom_types[parent.type]
            pterm = ff.polarizable_terms.get(atype.eqt_polar)
            if pterm is None:
                raise Exception(f'Polarizable term for {atype} not found in FF')
            if type(pterm) != DrudeTerm:
                raise Exception('Polarizable terms other than DrudeTerm haven\'t been implemented')
            drude.mass = pterm.mass
            parent.mass -= pterm.mass

    def assign_charge_from_ff(self, ff: ForceField):
        '''
        Assign charges of all atoms from the force field.
        The atom types should have been defined.
        The charge of corresponding AtomType in FF will be assigned to the atoms first.
        Then if there are ChargeIncrementTerm, the charge will be transferred between bonded atoms.
        The charge of Drude particles will be determined ONLY from the DrudePolarizableTerm,
        and the the same amount of charge will be subtracted from the parent atoms.
        That is, if there is an AtomType for Drude particles, the charge attribute of this AtomType will be simply ignored.
        '''
        _assigned = [False] * self.n_atom
        for atom in self._atoms:
            if atom.is_drude:
                continue
            try:
                charge = ff.atom_types[atom.type].charge
            except:
                raise Exception(f'Atom type {atom.type} not found in FF')
            atom.charge = charge
            if charge != 0:
                _assigned[atom.id_in_molecule] = True

        if len(ff.charge_increment_terms) > 0:
            for bond in filter(lambda x: not x.is_drude, self._bonds):
                try:
                    increment = ff.get_charge_increment(bond)
                except FFTermNotFoundError:
                    increment = 0
                bond.atom1.charge += increment
                bond.atom2.charge -= increment
                if increment != 0:
                    _assigned[bond.atom1.id_in_molecule] = True
                    _assigned[bond.atom2.id_in_molecule] = True

        for i, atom in enumerate(self.atoms):
            if not atom.is_drude and not _assigned[i]:
                logger.warning(f'Charge not assigned for {str(atom)} '
                               f'because both charge and increment in the FF are zero or not found')

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
            # store _alpha and _thole for PSF exporting
            parent._alpha = alpha
            parent._thole = pterm.thole
