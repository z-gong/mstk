import math
import copy
import numpy as np
from .atom import Atom
from .connectivity import *
from ..forcefield import FFSet
from ..forcefield.ffterm import *


class Molecule():
    def __init__(self, name: str = 'UNK'):
        self.id: int = -1  # id in topology
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

    def generate_angle_dihedral_improper(self):
        '''
        generate angle, dihedral and improper from bonds
        the existent angles, dihedrals and impropers will be replaced
        '''
        self._angles = []
        self._dihedrals = []
        self._impropers = []

        for bond in self._bonds:
            atom2 = bond.atom1
            atom3 = bond.atom2
            for atom1 in atom2.bond_partners:
                if atom1 != atom3:
                    self.add_angle(atom1, atom2, atom3, check_existence=True)
            for atom4 in atom3.bond_partners:
                if atom2 != atom4:
                    self.add_angle(atom2, atom3, atom4, check_existence=True)

            for atom1 in atom2.bond_partners:
                for atom4 in atom3.bond_partners:
                    if atom1 != atom3 and atom2 != atom4 and atom1 != atom4:
                        self.add_dihedral(atom1, atom2, atom3, atom4, check_existence=True)

        for atom in self._atoms:
            if len(atom.bond_partners) == 3:
                self.add_improper(atom, *atom.bond_partners)

    def guess_connectivity_from_ff(self, params: FFSet, bond_tolerance=0.025,
                                   angle_tolerance=15, pbc='xyz', box=None):
        '''
        Guess the connectivity from force field parameter set
        It requires that atoms types are defined and positions are provided
        The distance between nearby atoms will be calculated and compared with length in force field
        If the deviation is smaller than bond tolerance, the bond will be added
        Then the angles, dihedrals and impropers will be constructed from bonds
        All angles should be presented in force field parameter set
        If the deviation between angle and force field value is smaller than angle tolerance, the angle will be added
        As long as parameters for dihedrals and impropers present in force field, the dihedrals and impropers will be added
        '''
        if not self.has_position:
            raise Exception('Positions are required for guessing connectivity')

        self._bonds = []
        self._angles = []
        self._dihedrals = []
        self._impropers = []

        def check_and_add_angle(atom1, atom2, atom3):
            atype1 = params.atom_types[atom1.type].eqt_angle
            atype2 = params.atom_types[atom2.type].eqt_angle
            atype3 = params.atom_types[atom3.type].eqt_angle
            angle_term = AngleTerm(atype1, atype2, atype3, 0)
            if angle_term.name not in params.angle_terms.keys():
                raise Exception(f'{str(angle_term)} not exist in force field parameter set')
            angle_term: AngleTerm = params.angle_terms[angle_term.name]

            delta21 = atom1.position - atom2.position
            delta23 = atom3.position - atom2.position
            if 'x' in pbc:
                delta21[0] -= math.ceil(delta21[0] / box[0] - 0.5) * box[0]
                delta23[0] -= math.ceil(delta23[0] / box[0] - 0.5) * box[0]
            if 'y' in pbc:
                delta21[1] -= math.ceil(delta21[1] / box[1] - 0.5) * box[1]
                delta23[1] -= math.ceil(delta23[1] / box[1] - 0.5) * box[1]
            if 'z' in pbc:
                delta21[2] -= math.ceil(delta21[2] / box[2] - 0.5) * box[2]
                delta23[2] -= math.ceil(delta23[2] / box[2] - 0.5) * box[2]
            cos = delta21.dot(delta23) / delta21.dot(delta21) / delta23.dot(delta23)
            theta = np.arccos(np.clip(cos, -1, 1))
            if abs(theta * 180 / math.pi - angle_term.theta <= angle_tolerance):
                self.add_angle(atom1, atom2, atom3)
            else:
                print(f'warning: {str(Angle(atom1, atom2, atom3))} not added '
                      f'because far from equilibrium value')

        for i in range(self.n_atom):
            for j in range(i, self.n_atom):
                atom1 = self.atoms[i]
                atom2 = self.atoms[j]
                if atom1.type not in params.atom_types.keys():
                    raise Exception(
                        f'Atom type {atom1.type} not exist in force field parameter set')
                if atom2.type not in params.atom_types.keys():
                    raise Exception(
                        f'Atom type {atom2.type} not exist in force field parameter set')
                atype1 = params.atom_types[atom1.type].eqt_bond
                atype2 = params.atom_types[atom2.type].eqt_bond
                bond_term = BondTerm(atype1, atype2, 0)
                if bond_term.name not in params.bond_terms.keys():
                    raise Exception(f'{str(bond_term)} not exist in force field parameter set')
                bond_term: BondTerm = params.bond_terms[bond_term.name]

                delta = atom2.position - atom1.position
                if 'x' in pbc:
                    delta[0] -= math.ceil(delta[0] / box[0] - 0.5) * box[0]
                if 'y' in pbc:
                    delta[1] -= math.ceil(delta[1] / box[1] - 0.5) * box[1]
                if 'z' in pbc:
                    delta[2] -= math.ceil(delta[2] / box[2] - 0.5) * box[2]
                if any(delta - bond_term.length > bond_tolerance):
                    continue

                if abs(np.sqrt(delta.dot(delta)) - bond_term.length) <= bond_tolerance:
                    self.add_bond(atom1, atom2)

        for bond in self._bonds:
            atom2 = bond.atom1
            atom3 = bond.atom2
            for atom1 in atom2.bond_partners:
                if Angle(atom1, atom2, atom3) not in self._angles:
                    check_and_add_angle(atom1, atom2, atom3)
            for atom4 in atom3.bond_partners:
                if Angle(atom2, atom3, atom4) not in self._angles:
                    check_and_add_angle(atom2, atom3, atom4)

            for atom1 in atom2.bond_partners:
                for atom4 in atom3.bond_partners:
                    if atom1 == atom4 or Dihedral(atom1, atom2, atom3, atom4) in self._dihedrals:
                        continue
                    atype1 = params.atom_types[atom1.type].eqt_dihedral_side
                    atype2 = params.atom_types[atom2.type].eqt_dihedral_center
                    atype3 = params.atom_types[atom3.type].eqt_dihedral_center
                    atype4 = params.atom_types[atom4.type].eqt_dihedral_side
                    dihedral_term = DihedralTerm(atype1, atype2, atype3, atype4)
                    if dihedral_term.name not in params.dihedral_terms.keys():
                        print(f'warning: {str(Dihedral(atom1, atom2, atom3, atom4))} '
                              f'not added because not exist in force field parameter set')
                    else:
                        self.add_dihedral(atom1, atom2, atom3, atom4)

        for atom1 in self._atoms:
            if len(atom1.bond_partners) != 3:
                continue
            atom2, atom3, atom4 = atom1.bond_partners
            atype1 = params.atom_types[atom1.type].eqt_improper_center
            atype2 = params.atom_types[atom2.type].eqt_improper_side
            atype3 = params.atom_types[atom3.type].eqt_improper_side
            atype4 = params.atom_types[atom4.type].eqt_improper_side
            improper_term = ImproperTerm(atype1, atype2, atype3, atype4)
            if improper_term.name not in params.improper_terms.keys():
                print(f'warning: {str(Improper(atom1, atom2, atom3, atom4))} '
                      f'not added because not exist in force field parameter set')
            else:
                self.add_improper(atom1, atom2, atom3, atom4)

    def assign_charge_from_forcefield(self, params: FFSet):
        for atom in self._atoms:
            try:
                atom.charge = params.atom_types[atom.type].charge
            except:
                raise Exception(f'Atom type {atom.type} not exist in parameter set')

        for bond in self._bonds:
            atom1, atom2 = bond.atom1, bond.atom2
            try:
                atype1: AtomType = params.atom_types[atom1.type].eqt_charge_increment
                atype2: AtomType = params.atom_types[atom2.type].eqt_charge_increment
            except:
                raise Exception(f'Atom type {bond.atom1.type} or {bond.atom2.type}'
                                f'not exist in parameter set')
            binc = params.charge_increment_terms.get(
                ChargeIncrementTerm(atype1.name, atype2.name, 0).name)
            if binc is not None:
                direction = 1 if atype1.name == binc.type1 else -1
                atom1.charge += binc.value * direction
                atom2.charge -= binc.value * direction
