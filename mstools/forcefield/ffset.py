import math
from .ffterm import *


class FFSet():
    VDW_LONGRANGE_CORRECT = 'correct'
    VDW_LONGRANGE_SHIFT = 'shift'

    LJ_MIXING_NONE = 'none'
    LJ_MIXING_LB = 'lorentz-berthelot'
    LJ_MIXING_GEOMETRIC = 'geometric'

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

        self.vdw_cutoff = 1.2 # nm
        self.vdw_long_range = FFSet.VDW_LONGRANGE_CORRECT
        self.lj_mixing_rule = FFSet.LJ_MIXING_NONE
        self.scale_14_vdw = 1.0
        self.scale_14_coulomb = 1.0

    def get_vdw_term(self, type1: AtomType, type2: AtomType):
        '''
        Get vdW term between two atom types
        If not exist and it's LJ126 form, then generate it using combination rule
        '''
        at1 = type1.eqt_vdw
        at2 = type2.eqt_vdw
        vdw = VdwTerm(at1, at2)
        if at1 == at2:
            vdw = self.vdw_terms.get(vdw.name)
        else:
            vdw = self.pairwise_vdw_terms.get(vdw.name)

        if vdw is not None:
            return vdw

        lj1 = self.vdw_terms.get(VdwTerm(at1, at1).name)
        lj2 = self.vdw_terms.get(VdwTerm(at2, at2).name)
        if lj1 is None or lj2 is None:
            raise Exception(f'VdwTerm for {str(at1)} or {str(at2)} not found')
        if not isinstance(lj1, LJ126Term) or not isinstance(lj2, LJ126Term):
            raise Exception(f'Cannot using mixing rule for VdwTerm other than LJ126')

        if self.lj_mixing_rule == self.LJ_MIXING_LB:
            epsilon = math.sqrt(lj1.epsilon * lj2.epsilon)
            sigma = (lj1.sigma + lj2.sigma) / 2
        elif self.lj_mixing_rule == self.LJ_MIXING_GEOMETRIC:
            epsilon = math.sqrt(lj1.epsilon * lj2.epsilon)
            sigma = math.sqrt(lj1.sigma * lj2.sigma)
        else:
            raise Exception('Unknown mixing rule for LJ126 parameters')

        return LJ126Term(at1, at2, epsilon, sigma)

    @staticmethod
    def save_to(top, file):
        '''
        This method should be implemented by subclasses
        '''
        raise NotImplementedError('This method haven\'t been implemented')
