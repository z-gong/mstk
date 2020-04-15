import math
from ..ffterm import *


class FFSet():
    LJ_MIXING_NONE = 0
    LJ_MIXING_LB = 1
    LJ_MIXING_GEOMETRIC = 2

    def __init__(self):
        self.atom_types: {str: AtomType} = {}
        self.charge_increment_terms: {str: ChargeIncrementTerm} = {}
        self.bond_terms: {str: BondTerm} = {}
        self.angle_terms: {str: AngleTerm} = {}
        self.dihedral_terms: {str: DihedralTerm} = {}
        self.improper_terms: {str: ImproperTerm} = {}
        self.vdw_terms: {str: VdwTerm} = {}
        self.pairwise_vdw_terms: {str: VdwTerm} = {}
        self.polarizable_terms: {str: PolarizableTerm} = {}
        self.lj_mixing_rule = FFSet.LJ_MIXING_NONE
        self.scale_14_vdw = 1.0
        self.scale_14_coulomb = 1.0

    @property
    def use_only_lj126_for_vdw(self):
        '''
        Some logic can be simplified if only LJ126 is used for vdW interactions
        '''
        for vdw in tuple(self.vdw_terms.values()) + tuple(self.pairwise_vdw_terms.values()):
            if type(vdw) != LJ126Term:
                return False
        return True

    def generate_pairwise_lj126_terms(self):
        '''
        Explicitly generate pairwise LJ126 terms using combination rule
        The ones already exist will not be touched
        '''
        atom_types = list(self.atom_types.values())
        for i in range(len(atom_types)):
            for j in range(i + 1, len(atom_types)):
                vdw = self.get_vdw_term(atom_types[i], atom_types[j])
                self.vdw_terms[vdw.name] = vdw

    def get_vdw_term(self, type1: AtomType, type2: AtomType):
        '''
        Get vdW term between two atom types
        If not exist and it's LJ126 form, then generate it using combination rule
        '''
        atype1 = type1.eqt_vdw
        atype2 = type2.eqt_vdw
        if atype1 == atype2:
            vdw = self.vdw_terms.get(VdwTerm(atype1.name, atype2.name).name)
        else:
            vdw = self.pairwise_vdw_terms.get(VdwTerm(atype1.name, atype2.name).name)

        if vdw is not None:
            return vdw

        lj1 = self.vdw_terms.get(VdwTerm(atype1.name, atype1.name).name)
        lj2 = self.vdw_terms.get(VdwTerm(atype2.name, atype2.name).name)
        if lj1 is None or lj2 is None:
            raise Exception(f'VdwTerm for {str(atype1)} or {str(atype2)} not found')
        if not isinstance(lj1, LJ126Term) or not isinstance(lj2, LJ126Term):
            raise Exception(f'Cannot using combination rule for VdwTerm other than LJ126')

        if self.lj_mixing_rule == self.LJ_MIXING_LB:
            epsilon = math.sqrt(lj1.epsilon * lj2.epsilon)
            sigma = (lj1.sigma + lj2.sigma) / 2
        elif self.lj_mixing_rule == self.LJ_MIXING_GEOMETRIC:
            epsilon = math.sqrt(lj1.epsilon * lj2.epsilon)
            sigma = math.sqrt(lj1.sigma * lj2.sigma)
        else:
            raise Exception('Unknown mixing rule for LJ126 parameters')

        return LJ126Term(atype1.name, atype2.name, epsilon, sigma)
