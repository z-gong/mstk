import math
from typing import Union
from .ffterm import *
from .errors import *


class ForceField():
    '''
    ForceField is a set of FFTerms describing the interactions between atoms in a topology.

    Attributes
    ----------
    atom_types : dict, [str, AtomType]

    '''
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

        self.vdw_cutoff = 1.2  # nm
        self.vdw_long_range = ForceField.VDW_LONGRANGE_CORRECT
        self.lj_mixing_rule = ForceField.LJ_MIXING_NONE
        self.scale_14_vdw = 1.0
        self.scale_14_coulomb = 1.0

    @staticmethod
    def open(*files):
        '''
        Load ForceField from files.

        The parser for reading the files will be determined from the extension of the first file.
        Therefore, all the files should be in the same format.

        Parameters
        ----------
        files : list of str
            The files to be read. The extension can be `zfp`, `ppf` or `ff`.
            `zfp` file is the native format of mstools, and will be parsed by Zfp.
            `ppf` file will be parsed as DFF format by Ppf.
            `ff` file will be parsed as fftool format by Padua.

        Returns
        -------
        forcefield : ForceField

        '''
        from .zfp import Zfp
        from .ppf import Ppf
        from .padua import Padua

        file = files[0]
        if file.endswith('.ff'):
            return Padua(*files).forcefield
        elif file.endswith('.ppf'):
            return Ppf(*files).forcefield
        elif file.endswith('.zfp'):
            return Zfp(*files).forcefield
        else:
            raise Exception('Unsupported format')

    def write(self, file):
        from .zfp import Zfp
        from .ppf import Ppf
        from .padua import Padua

        if file.endswith('.zfp'):
            Zfp.save_to(self, file)
        elif file.endswith('.ppf'):
            Ppf.save_to(self, file)
        else:
            raise Exception('Unsupported format for forcefield output')

    @staticmethod
    def save_to(ff, file):
        '''
        This method should be implemented by subclasses
        '''
        raise NotImplementedError('This method haven\'t been implemented')

    @property
    def is_polarizable(self) -> bool:
        return len(self.polarizable_terms) != 0

    def get_settings(self) -> {str: str}:
        d = {}
        for i in ('vdw_cutoff', 'vdw_long_range', 'lj_mixing_rule',
                  'scale_14_vdw', 'scale_14_coulomb'):
            d[i] = getattr(self, i)
        return d

    def restore_settings(self, d):
        for i in ('vdw_cutoff', 'vdw_long_range', 'lj_mixing_rule',
                  'scale_14_vdw', 'scale_14_coulomb'):
            setattr(self, i, d[i])

    def add_term(self, term: FFTerm, replace=False):
        _duplicated = False
        if isinstance(term, AtomType):
            if term.name not in self.atom_types or replace:
                self.atom_types[term.name] = term
            else:
                _duplicated = True
        elif isinstance(term, BondTerm):
            if term.name not in self.bond_terms or replace:
                self.bond_terms[term.name] = term
            else:
                _duplicated = True
        elif isinstance(term, AngleTerm):
            if term.name not in self.angle_terms or replace:
                self.angle_terms[term.name] = term
            else:
                _duplicated = True
        elif isinstance(term, DihedralTerm):
            if term.name not in self.dihedral_terms or replace:
                self.dihedral_terms[term.name] = term
            else:
                _duplicated = True
        elif isinstance(term, ImproperTerm):
            if term.name not in self.improper_terms or replace:
                self.improper_terms[term.name] = term
            else:
                _duplicated = True
        elif isinstance(term, VdwTerm):
            if term.type1 == term.type2:
                if term.name not in self.vdw_terms or replace:
                    self.vdw_terms[term.name] = term
                else:
                    _duplicated = True
            else:
                if term.name not in self.pairwise_vdw_terms or replace:
                    self.pairwise_vdw_terms[term.name] = term
                else:
                    _duplicated = True
        elif isinstance(term, PolarizableTerm):
            if term.name not in self.polarizable_terms or replace:
                self.polarizable_terms[term.name] = term
            else:
                _duplicated = True
        else:
            raise Exception('Invalid term to add')

        if _duplicated:
            raise Exception(f'{str(term)} already exist in FF')

    def get_vdw_term(self, type1: Union[AtomType, str], type2: Union[AtomType, str],
                     mixing=True) -> VdwTerm:
        '''
        Get vdW term between two atom types
        If not exist and mixing=True and it's LJ126 form, then generate it using mixing rule
        '''
        at1: str = type1.eqt_vdw if type(type1) is AtomType else type1
        at2: str = type2.eqt_vdw if type(type2) is AtomType else type2
        vdw = VdwTerm(at1, at2)
        if at1 == at2:
            vdw = self.vdw_terms.get(vdw.name)
            if vdw is not None:
                return vdw
            else:
                raise FFTermNotFoundError(f'VdwTerm for {at1} not found')

        vdw = self.pairwise_vdw_terms.get(vdw.name)

        if vdw is not None:
            return vdw
        elif not mixing:
            raise FFTermNotFoundError(f'VdwTerm for {at1}-{at2} not found')

        lj1 = self.vdw_terms.get(VdwTerm(at1, at1).name)
        lj2 = self.vdw_terms.get(VdwTerm(at2, at2).name)
        if lj1 is None or lj2 is None:
            raise FFTermNotFoundError(f'VdwTerm for {at1} or {at2} not found')
        if type(lj1) is not LJ126Term or type(lj2) is not LJ126Term:
            raise Exception('Cannot use mixing rule for vdW terms other than LJ126Term')

        if self.lj_mixing_rule == self.LJ_MIXING_LB:
            epsilon = math.sqrt(lj1.epsilon * lj2.epsilon)
            sigma = (lj1.sigma + lj2.sigma) / 2
        elif self.lj_mixing_rule == self.LJ_MIXING_GEOMETRIC:
            epsilon = math.sqrt(lj1.epsilon * lj2.epsilon)
            sigma = math.sqrt(lj1.sigma * lj2.sigma)
        else:
            raise Exception('Unknown mixing rule for LJ126Term')

        return LJ126Term(at1, at2, epsilon, sigma)

    def get_eqt_for_bond(self, bond) -> (str, str):
        try:
            at1 = self.atom_types.get(bond.atom1.type).eqt_bond
            at2 = self.atom_types.get(bond.atom2.type).eqt_bond
        except:
            raise FFTermNotFoundError(
                f'Atom type {bond.atom1.type} or {bond.atom2.type} not found in FF')
        return at1, at2

    def get_bond_term(self, type1, type2) -> BondTerm:
        at1: str = type1.eqt_bond if type(type1) is AtomType else type1
        at2: str = type2.eqt_bond if type(type2) is AtomType else type2
        bterm = BondTerm(at1, at2, 0)
        try:
            bterm = self.bond_terms[bterm.name]
        except:
            raise FFTermNotFoundError(f'{str(bterm)} not found in FF')

        return bterm

    def get_eqt_for_angle(self, angle) -> (str, str, str):
        try:
            at1 = self.atom_types.get(angle.atom1.type).eqt_ang_s
            at2 = self.atom_types.get(angle.atom2.type).eqt_ang_c
            at3 = self.atom_types.get(angle.atom3.type).eqt_ang_s
        except:
            raise FFTermNotFoundError(f'Atom type {angle.atom1.type} or {angle.atom2.type} '
                                      f'or {angle.atom3.type} not found in FF')
        return at1, at2, at3

    def get_angle_term(self, type1, type2, type3) -> AngleTerm:
        at1: str = type1.eqt_ang_s if type(type1) is AtomType else type1
        at2: str = type2.eqt_ang_c if type(type2) is AtomType else type2
        at3: str = type3.eqt_ang_s if type(type3) is AtomType else type3
        aterm = AngleTerm(at1, at2, at3, 0)
        try:
            aterm = self.angle_terms[aterm.name]
        except:
            raise FFTermNotFoundError(f'{str(aterm)} not found in FF')

        return aterm

    def get_eqt_for_dihedral(self, dihedral) -> [(str, str, str, str)]:
        try:
            at1 = self.atom_types.get(dihedral.atom1.type).eqt_dih_s
            at2 = self.atom_types.get(dihedral.atom2.type).eqt_dih_c
            at3 = self.atom_types.get(dihedral.atom3.type).eqt_dih_c
            at4 = self.atom_types.get(dihedral.atom4.type).eqt_dih_s
        except:
            raise FFTermNotFoundError(f'Atom type {dihedral.atom1.type} or {dihedral.atom2.type} '
                                      f'or {dihedral.atom3.type} or {dihedral.atom4.type} not found in FF')
        return [(at1, at2, at3, at4),
                (at1, at2, at3, '*'),
                ('*', at2, at3, at4),
                ('*', at2, at3, '*')]

    def get_dihedral_term(self, type1, type2, type3, type4) -> DihedralTerm:
        at1: str = type1.eqt_dih_s if type(type1) is AtomType else type1
        at2: str = type2.eqt_dih_c if type(type2) is AtomType else type2
        at3: str = type3.eqt_dih_c if type(type3) is AtomType else type3
        at4: str = type4.eqt_dih_s if type(type4) is AtomType else type4
        dterm = DihedralTerm(at1, at2, at3, at4)
        try:
            dterm = self.dihedral_terms[dterm.name]
        except:
            raise FFTermNotFoundError(f'{str(dterm)} with or w/o wildcard not found in FF')

        return dterm

    def get_eqt_for_improper(self, improper) -> [(str, str, str, str)]:
        try:
            at1 = self.atom_types.get(improper.atom1.type).eqt_imp_c
            at2 = self.atom_types.get(improper.atom2.type).eqt_imp_s
            at3 = self.atom_types.get(improper.atom3.type).eqt_imp_s
            at4 = self.atom_types.get(improper.atom4.type).eqt_imp_s
        except:
            raise FFTermNotFoundError(f'Atom type {improper.atom1.type} or {improper.atom2.type} '
                                      f'or {improper.atom3.type} or {improper.atom4.type} not found in FF')
        return [(at1, at2, at3, at4),
                (at1, at2, at3, '*'),
                (at1, at2, '*', at4),
                (at1, '*', at3, at4),
                (at1, at2, '*', '*'),
                (at1, '*', at3, '*'),
                (at1, '*', '*', at4),
                (at1, '*', '*', '*')]

    def get_improper_term(self, type1, type2, type3, type4) -> ImproperTerm:
        at1: str = type1.eqt_imp_c if type(type1) is AtomType else type1
        at2: str = type2.eqt_imp_s if type(type2) is AtomType else type2
        at3: str = type3.eqt_imp_s if type(type3) is AtomType else type3
        at4: str = type4.eqt_imp_s if type(type4) is AtomType else type4
        iterm = ImproperTerm(at1, at2, at3, at4)
        try:
            iterm = self.improper_terms[iterm.name]
        except:
            raise FFTermNotFoundError(f'{str(iterm)} with or w/o wildcard not found in FF')

        return iterm

    def get_charge_increment(self, bond) -> float:
        try:
            at1 = self.atom_types.get(bond.atom1.type).eqt_q_inc
            at2 = self.atom_types.get(bond.atom2.type).eqt_q_inc
        except:
            raise FFTermNotFoundError(
                f'AtomType {bond.atom1.type} or {bond.atom2.type} not found in FF')
        qterm = ChargeIncrementTerm(at1, at2, 0)
        try:
            qterm = self.charge_increment_terms[qterm.name]
        except:
            raise FFTermNotFoundError(f'{str(qterm)} not found in FF')

        direction = 1 if at1 == qterm.type1 else -1
        return qterm.value * direction
