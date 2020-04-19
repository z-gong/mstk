import math
import random
import copy
import warnings
import numpy as np
from .atom import Atom
from .connectivity import *
from ..forcefield import FFSet
from ..forcefield.ffterm import *


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

    def __str__(self):
        return f'<Molecule: {self.name} {self.id}>'

    def __repr__(self):
        return str(self) + ' at 0x' + str(hex(id(self))[2:].upper())

    def __deepcopy__(self, memodict={}):
        mol = Molecule()
        mol.id = self.id
        mol.name = self.name
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

    @property
    def topology(self):
        return self._topology

    def add_atom(self, atom: Atom):
        self._atoms.append(atom)
        atom._molecule = self
        atom._id_in_molecule = len(self._atoms) - 1
        if self._topology is not None:
            self._topology.init_from_molecules(self._topology.molecules, deepcopy=False)

    def remove_atom(self, atom: Atom):
        self._atoms.remove(atom)
        atom._molecule = None
        atom._id_in_molecule = -1
        if self._topology is not None:
            self._topology.init_from_molecules(self._topology.molecules, deepcopy=False)
        for bond in atom._bonds[:]:
            self.remove_bond(bond)

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
    def initialized(self):
        return self.id != -1

    @property
    def has_position(self):
        return all(atom.has_position for atom in self.atoms)

    def get_drude_pairs(self)->[(Atom, Atom)]:
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

    def generate_angle_dihedral_improper(self):
        '''
        Generate angle, dihedral and improper from bonds
        The existing angles, dihedrals and impropers will be replaced
        The Drude particles will be ignored
        '''
        self._angles = []
        self._dihedrals = []
        self._impropers = []

        for bond in filter(lambda x: not x.is_drude, self._bonds):
            atom2 = bond.atom1
            atom3 = bond.atom2
            for atom1 in filter(lambda x: not x.is_drude, atom2.bond_partners):
                if atom1 != atom3:
                    self.add_angle(atom1, atom2, atom3, check_existence=True)
            for atom4 in filter(lambda x: not x.is_drude, atom3.bond_partners):
                if atom2 != atom4:
                    self.add_angle(atom2, atom3, atom4, check_existence=True)

            for atom1 in filter(lambda x: not x.is_drude, atom2.bond_partners):
                for atom4 in filter(lambda x: not x.is_drude, atom3.bond_partners):
                    if atom1 != atom3 and atom2 != atom4 and atom1 != atom4:
                        self.add_dihedral(atom1, atom2, atom3, atom4, check_existence=True)

        for atom in self._atoms:
            if len([p for p in atom.bond_partners if not p.is_drude]) == 3:
                self.add_improper(atom, *atom.bond_partners)

    def guess_connectivity_from_forcefield(self, params: FFSet, bond_length_limit=0.25,
                                           bond_tolerance=0.025, angle_tolerance=None,
                                           pbc='xyz', box=None):
        '''
        Guess bonds, angles, dihedrals and impropers from force field.
        It requires that atoms types are defined and positions are available.
        The distance between nearby atoms will be calculated.
        If it's smaller than bond_length_limit, then it will be compared with the equilibrium length in FF.
        If the deviation is smaller than bond_tolerance, the bond will be added.
        Then angles will be constructed from bonds. If angle_tolerance is None, all angles will be added.
        If angle_tolerance is set (as degree), then the angles will be compared with the equilibrium value in FF,
        and be added only if the deviation is smaller than angle_tolerance.
        Dihedrals and impropers will be constructed form bonds and be added if relevant terms are presented in FF.
        '''
        if not self.has_position:
            raise Exception('Positions are required for guessing connectivity')

        if any(atom.is_drude for atom in self._atoms):
            raise Exception('Drude particles should be removed before guess connectivity')

        self._bonds = []
        self._angles = []
        self._dihedrals = []
        self._impropers = []

        for i in range(self.n_atom):
            for j in range(i, self.n_atom):
                atom1 = self.atoms[i]
                atom2 = self.atoms[j]
                try:
                    at1 = params.atom_types[atom1.type].eqt_bond
                    at2 = params.atom_types[atom2.type].eqt_bond
                except:
                    raise Exception(f'AtomType {atom1.type} or {atom2.type} not found in FF')

                delta = atom2.position - atom1.position
                if 'x' in pbc and box is not None:
                    delta[0] -= math.ceil(delta[0] / box[0] - 0.5) * box[0]
                if 'y' in pbc and box is not None:
                    delta[1] -= math.ceil(delta[1] / box[1] - 0.5) * box[1]
                if 'z' in pbc and box is not None:
                    delta[2] -= math.ceil(delta[2] / box[2] - 0.5) * box[2]
                if any(delta > bond_length_limit):
                    continue

                bterm = BondTerm(at1, at2, 0)
                if bterm.name not in params.bond_terms.keys():
                    warnings.warn(
                        f'Distance between {str(atom1)} and {str(atom2)} smaller than limit '
                        f'but {str(bterm)} not found in FF')
                    continue

                bterm: BondTerm = params.bond_terms[bterm.name]
                if any(delta - bterm.length > bond_tolerance):
                    continue

                if abs(np.sqrt(delta.dot(delta)) - bterm.length) <= bond_tolerance:
                    self.add_bond(atom1, atom2)

        # generate angles etc..., and then remove them if requirements are not satisfied
        self.generate_angle_dihedral_improper()

        if angle_tolerance is not None:
            for angle in self._angles[:]:
                at1 = params.atom_types[angle.atom1.type].eqt_ang_s
                at2 = params.atom_types[angle.atom2.type].eqt_ang_c
                at3 = params.atom_types[angle.atom3.type].eqt_ang_s
                aterm = AngleTerm(at1, at2, at3, 0)
                if aterm.name not in params.angle_terms.keys():
                    raise Exception(f'{str(angle)} constructed but {str(aterm)} not found in FF')

                aterm: AngleTerm = params.angle_terms[aterm.name]
                delta21 = angle.atom1.position - angle.atom2.position
                delta23 = angle.atom3.position - angle.atom2.position
                if 'x' in pbc and box is not None:
                    delta21[0] -= math.ceil(delta21[0] / box[0] - 0.5) * box[0]
                    delta23[0] -= math.ceil(delta23[0] / box[0] - 0.5) * box[0]
                if 'y' in pbc and box is not None:
                    delta21[1] -= math.ceil(delta21[1] / box[1] - 0.5) * box[1]
                    delta23[1] -= math.ceil(delta23[1] / box[1] - 0.5) * box[1]
                if 'z' in pbc and box is not None:
                    delta21[2] -= math.ceil(delta21[2] / box[2] - 0.5) * box[2]
                    delta23[2] -= math.ceil(delta23[2] / box[2] - 0.5) * box[2]
                cos = delta21.dot(delta23) / delta21.dot(delta21) / delta23.dot(delta23)
                theta = np.arccos(np.clip(cos, -1, 1))
                if abs(theta * 180 / math.pi - aterm.theta) > angle_tolerance:
                    self.remove_angle(angle)
                    warnings.warn(f'{str(angle)} not added because '
                                  f'its value {theta} is far from equilibrium {aterm.theta}')

        for dihedral in self._dihedrals[:]:
            at1 = params.atom_types[dihedral.atom1.type].eqt_dih_s
            at2 = params.atom_types[dihedral.atom2.type].eqt_dih_c
            at3 = params.atom_types[dihedral.atom3.type].eqt_dih_c
            at4 = params.atom_types[dihedral.atom4.type].eqt_dih_s
            dterm = DihedralTerm(at1, at2, at3, at4)
            if dterm.name not in params.dihedral_terms.keys():
                self.remove_dihedral(dihedral)
                warnings.warn(f'{str(dihedral)} not added because {str(dterm)} not found in FF')

        for improper in self._impropers[:]:
            at1 = params.atom_types[improper.atom1.type].eqt_imp_c
            at2 = params.atom_types[improper.atom2.type].eqt_imp_s
            at3 = params.atom_types[improper.atom3.type].eqt_imp_s
            at4 = params.atom_types[improper.atom4.type].eqt_imp_s
            iterm = ImproperTerm(at1, at2, at3, at4)
            if iterm.name not in params.improper_terms.keys():
                self.remove_improper(improper)
                warnings.warn(f'{str(improper)} not added because {str(iterm)} not foundn in FF')

    def generate_drude_particles(self, params: FFSet, type_drude='DRUDE', seed=1):
        '''
        Generate Drude particles from polarizable terms in force field.
        The atom types should have been assigned already.
        Drude particle will not be generated if polarizable term for its atom type can not be found in the FF.
        Note that The existing Drude particles will be deleted before generating.
        The mass defined in the FF term will be transferred from parent atom to the Drude particle.
        The Drude charge will be calculated from the FF term and transferred from parent atom to the Drude particle.
        Bonds between parent-Drude will be generated and added to the topology.
        '''
        if len(params.polarizable_terms) == 0:
            raise Exception('PolarizableTerms not found in force field')

        random.seed(seed)

        self.remove_drude_particles()
        for atom in self._atoms[:]:
            atype = params.atom_types.get(atom.type)
            if atype is None:
                raise Exception(f'Atom type {atom.type} not found in force field')

            pterm = params.polarizable_terms.get(atype.eqt_polar)
            if pterm is None:
                continue
            if type(pterm) != IsotropicDrudeTerm:
                raise Exception('Polarizable terms other than IsotropicDrudeTerm '
                                'haven\'t been implemented')
            drude = Atom()
            drude.is_drude = True
            drude.type = type_drude
            drude.mass = pterm.mass
            atom.mass -= drude.mass
            drude.charge = - pterm.get_charge()
            atom.charge += pterm.get_charge()
            if atom.has_position:
                deviation = random.random() / 1000
                drude.position = atom.position + [deviation, deviation, deviation]
                drude.has_position = True
            self._atoms.append(drude)
            self.add_bond(atom, drude)
        if self._topology is not None:
            self._topology.init_from_molecules(self._topology.molecules, deepcopy=False)

        dtype = params.atom_types.get(type_drude)
        if dtype is None:
            dtype = AtomType(type_drude)
            params.atom_types[dtype.name] = dtype
            warnings.warn(f'AtomType for Drude particle not found in FF. '
                          f'{str(dtype)} is added to the FF')
        vdw = params.vdw_terms.get(dtype.eqt_vdw)
        if vdw is None:
            vdw = LJ126Term(dtype.eqt_vdw, dtype.eqt_vdw, 0.0, 1.0)
            params.vdw_terms[vdw.name] = vdw
            warnings.warn(f'VdwTerm for Drude particle not found in FF. '
                          f'{str(vdw)} with zero interactions is added to the FF')

    def remove_drude_particles(self):
        for atom in self._atoms[:]:
            if atom.is_drude:
                self.remove_bond(atom._bonds[0])
                self._atoms.remove(atom)
        if self._topology is not None:
            self._topology.init_from_molecules(self._topology.molecules, deepcopy=False)

    def assign_mass_from_forcefield(self, params: FFSet):
        '''
        Assign masses of all atoms from the force field.
        The atom types should have been assigned already, and the AtomType in FF should carry mass information.
        The masses of Drude particles will be determined ONLY from the IsotropicDrudeTerm,
        and the the same amount of mass will be subtracted from the parent atoms.
        That is, if there is an AtomType for Drude particles, the mass attribute of this AtomType will be simply ignored.
        '''
        for atom in self._atoms:
            if atom.is_drude:
                continue
            atype = params.atom_types.get(atom.type)
            if atype is None:
                raise Exception(f'Atom type {atom.type} not found in force field')
            if atype.mass == -1:
                raise Exception(f'{atype} does not carry mass information')
            if atype.mass == 0:
                warnings.warn(f'{atype} carries zero mass. You should make sure it is correct')
            atom.mass = atype.mass

        for parent, drude in self.get_drude_pairs():
            atype = params.atom_types[parent.type]
            pterm = params.polarizable_terms.get(atype.eqt_polar)
            if pterm is None:
                raise Exception(f'Polarizable term for {atype} not found in FF')
            if type(pterm) != IsotropicDrudeTerm:
                raise Exception('Polarizable terms other than IsotropicDrudeTerm '
                                'haven\'t been implemented')
            drude.mass = pterm.mass
            parent.mass -= pterm.mass

    def assign_charge_from_forcefield(self, params: FFSet):
        '''
        Assign charges of all atoms from the force field.
        The atom types should have been assigned already.
        The charge of corresponding AtomType in FF will be assigned to the atoms first.
        Then if there are ChargeIncrementTerm, the charge will be transferred between bonded atoms.
        The charge of Drude particles will be determined ONLY from the IsotropicDrudeTerm,
        and the the same amount of charge will be subtracted from the parent atoms.
        That is, if there is an AtomType for Drude particles, the charge attribute of this AtomType will be simply ignored.
        '''
        for atom in self._atoms:
            if atom.is_drude:
                continue
            try:
                atom.charge = params.atom_types[atom.type].charge
            except:
                raise Exception(f'Atom type {atom.type} not found in force field')

        if len(params.charge_increment_terms) > 0:
            for bond in self._bonds:
                if bond.is_drude:
                    continue
                atom1, atom2 = bond.atom1, bond.atom2
                at1 = params.atom_types[atom1.type].eqt_q_inc
                at2 = params.atom_types[atom2.type].eqt_q_inc
                binc = ChargeIncrementTerm(at1, at2, 0)
                try:
                    binc = params.charge_increment_terms[binc.name]
                except:
                    warnings.warn(f'{binc} not found in force field')

                direction = 1 if at1.name == binc.type1 else -1
                atom1.charge += binc.value * direction
                atom2.charge -= binc.value * direction

        for parent, drude in self.get_drude_pairs():
            atype = params.atom_types[parent.type]
            pterm = params.polarizable_terms.get(atype.eqt_polar)
            if pterm is None:
                raise Exception(f'Polarizable term for {atype} not found in FF')
            if type(pterm) != IsotropicDrudeTerm:
                raise Exception('Polarizable terms other than IsotropicDrudeTerm '
                                'haven\'t been implemented')
            drude.charge = -pterm.get_charge()
            parent.charge += pterm.get_charge()
