import math
from typing import Union
from .ffterm import *
from .errors import *


class ForceField():
    '''
    ForceField is a set of :class:`FFTerm` objects describing the energy function of
    a :class:`~mstools.topology.Topology` with regard to the positions of atoms.

    Attributes
    ----------
    atom_types : dict, [str, AtomType]
        All the atom types in this force field.
        For each key value pair, the value is a AtomType object,
        the key is the :attr:`~AtomType.name` of the object.
    bci_terms : dict, [str, ChargeIncrementTerm]
        All the charge increment terms.
        For each key value pair, the value is a ChargeIncrementTerm object,
        the key is the :attr:`~ChargeIncrementTerm.name` of the object.
    vdw_terms : dict, [str, subclass of VdwTerm]
        vdW terms between same types.
        For each key value pair, the value is a object of subclass of VdwTerm.
        the key is the :attr:`~VdwTerm.name` of the object.
        The term can be either :class:`LJ126Term` or :class:`MieTerm`,
        as long as it's a subclass of VdwTerm.
        Because we are using name of the term as key,
        there can only be one term to describe the vdW interactions between two atom types.
        Which means, there cannot be a LJ126Term and a MieTerm with same name in a ForceField.
        However, you can have LJ126Term for one pair of atom types and MieTerm for another pairs.
    pairwise_vdw_terms : dict, [str, subclass of VdwTerm]
        Pairwise vdW terms between unlike types.
        For each key value pair, the value is a object of subclass of VdwTerm.
        the key is the :attr:`~VdwTerm.name` of the object.
        If a pair of unlike atom types is not included in this dict,
        then the vdW interactions between this pair will be generated from mixing rule.
    polarizable_terms : dict, [str, subclass of PolarizableTerm]
        Polarization terms for polarizable atom types.
        For each key value pair, the value is a object of subclass of PolarizableTerm.
        the key is the :attr:`~PolarizableTerm.name` of the object.
        Currently, isotropic Drude model is implemented.
    bond_terms : dict, [str, subclass of BondTerm]
        Bond terms between two atom types.
        For each key value pair, the value is a object of subclass of BondTerm.
        the key is the :attr:`~BondTerm.name` of the object.
    angle_terms : dict, [str, subclass of AngleTerm]
        Angle terms between three atom types.
        For each key value pair, the value is a object of subclass of AngleTerm.
        the key is the :attr:`~AngleTerm.name` of the object.
    dihedral_terms : dict, [str, subclass of DihedralTerm]
        Dihedral terms between four atom types.
        For each key value pair, the value is a object of subclass of DihedralTerm.
        the key is the :attr:`~DihedralTerm.name` of the object.
    improper_terms : dict, [str, subclass of ImproperTerm]
        Improper terms between four atom types.
        For each key value pair, the value is a object of subclass of ImproperTerm.
        the key is the :attr:`~ImproperTerm.name` of the object.
    vdw_cutoff : float
        Cutoff distance for vdW interactions.
        Will be ignored if periodic boundary condition is not used.
    vdw_long_range : str
        The method to handle the long-range vdW interactions.
        Optional values: :attr:`VDW_LONGRANGE_CORRECT`, :attr:`VDW_LONGRANGE_SHIFT`
        Will be ignored if periodic boundary condition is not used.
    lj_mixing_rule : str
        The rule to generate the pairwise LJ126 parameters if not provided.
        Optional values: :attr:`LJ_MIXING_NONE`, :attr:`LJ_MIXING_LB`, :attr:`LJ_MIXING_GEOMETRIC`
    scale_14_vdw : float
        The scaling factor for 1-4 vdW interactions
    scale_14_coulomb : float
        The scaling factor for 1-4 Coulomb interactions
    '''
    #: Truncated vdW interactions with long range continuum correction
    VDW_LONGRANGE_CORRECT = 'correct'
    #: Truncated vdW interactions with potential shifted so that vdW energy equal to zero at cutoff
    VDW_LONGRANGE_SHIFT = 'shift'

    #: Do not use mixing rule for LJ126 parameters
    LJ_MIXING_NONE = 'none'
    #: Lorentz-Berthelot mixing rule for LJ126 parameters (used by CHARMM and AMBER)
    LJ_MIXING_LB = 'lorentz-berthelot'
    #: Geometric mixing rule for LJ126 parameters (used by OPLS)
    LJ_MIXING_GEOMETRIC = 'geometric'

    def __init__(self):
        self.atom_types: {str: AtomType} = {}
        self.bci_terms: {str: ChargeIncrementTerm} = {}
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

            * `zfp` file is the native format of mstools, and will be parsed by :class:`Zfp`.
            * `ppf` file will be parsed as DFF format by :class:`Ppf`.
            * `ff` file will be parsed as fftool format by :class:`Padua`.

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
        '''
        Write this force field into a file.

        The format of the file will be determined by its extension

        Parameters
        ----------
        file : str
        '''

        from .zfp import Zfp
        from .ppf import Ppf
        from .padua import Padua

        if file.endswith('.zfp'):
            Zfp.save_to(self, file)
        elif file.endswith('.ppf'):
            Ppf.save_to(self, file)
        else:
            raise Exception('Unsupported format for forcefield output')

    @property
    def is_polarizable(self):
        '''
        Whether or not this is a polarizable force field.

        Return True as long as the :attr:`polarizable_terms` is not empty.

        Returns
        -------
        is : bool
        '''
        return len(self.polarizable_terms) != 0

    def get_settings(self):
        '''
        Get the settings of this force field.

        Returns
        -------
        settings : dict, [str, str]

        '''
        d = {}
        for i in ('vdw_cutoff', 'vdw_long_range', 'lj_mixing_rule',
                  'scale_14_vdw', 'scale_14_coulomb'):
            d[i] = getattr(self, i)
        return d

    def restore_settings(self, d):
        '''
        Restore the settings of this force field

        Parameters
        ----------
        d : dict, [str, str]

        '''
        for i in ('vdw_cutoff', 'vdw_long_range', 'lj_mixing_rule',
                  'scale_14_vdw', 'scale_14_coulomb'):
            setattr(self, i, d[i])

    def add_term(self, term, replace=False):
        '''
        Add a force field term into this force field.

        If a terms already exists in the force field,
        the argument `replace` will determine the behavior.
        If replace is True, the one existed will be replaced by the new term.
        Otherwise, an Exception will be raised.

        Parameters
        ----------
        term : subclass of FFTerm
        replace : bool
        '''
        _duplicated = False
        if isinstance(term, AtomType):
            if term.name not in self.atom_types or replace:
                self.atom_types[term.name] = term
            else:
                _duplicated = True
        elif isinstance(term, ChargeIncrementTerm):
            if term.name not in self.bci_terms or replace:
                self.bci_terms[term.name] = term
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

    def get_vdw_term(self, type1, type2, mixing=True):
        '''
        Get the vdW term between two atom types.

        If type1 and type2 are the same, then will search parameters in :attr:`vdw_terms` and return the reference of the term.
        If not found, an Exception will be raised.
        If type1 and type2 are not the same, then will search parameters in :attr:`pairwise_vdw_terms` and return the reference of the term.
        If not found and `mixing` is True and :attr:`lj_mixing_rule` does not equal to :attr:`LJ_MIXING_NONE`,
        and the vdw terms for both type1 and type2 are :class:`LJ126Term`,
        then a LJ126Term will be constructed from the mixing rule and be returned.
        Otherwise, an Exception will be raised.

        Parameters
        ----------
        type1 : AtomType or str
            If type1 is a AtomType object, then will search for its equivalent type for vdW parameters.
            If type1 is a str, then will match exactly this type.
        type2 : AtomType or str
            If type2 is a AtomType object, then will search for its equivalent type for vdW parameters.
            If type2 is a str, them will match exactly this type.
        mixing : bool
            Whether or not using combination rule for LJ126Term.

        Returns
        -------
        term : subclass of VdwTerm
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

    def get_eqt_for_bond(self, bond):
        '''
        Get the equivalent atom types for bond parameters of the two atoms in a bond.

        Parameters
        ----------
        bond : ~mstools.topology.Bond

        Returns
        -------
        type1 : str
            The equivalent atom type for bond parameters of the first atom in this bond.
        type2 : str
            The equivalent atom type for bond parameters of the second atom in this bond.
        '''
        try:
            at1 = self.atom_types.get(bond.atom1.type).eqt_bond
            at2 = self.atom_types.get(bond.atom2.type).eqt_bond
        except:
            raise FFTermNotFoundError(
                f'Atom type {bond.atom1.type} or {bond.atom2.type} not found in FF')
        return at1, at2

    def get_bond_term(self, type1, type2):
        '''
        Get the bond term for the bond formed by two atom types.

        If not found, an Exception will be raised.

        Parameters
        ----------
        type1 : AtomType or str
            If type1 is a AtomType object, then will search for its equivalent type for bond parameters.
            If type1 is a str, then will match exactly this type.
        type2 : AtomType or str
            If type2 is a AtomType object, then will search for its equivalent type for bond parameters.
            If type2 is a str, them will match exactly this type.

        Returns
        -------
        term : subclass of BondTerm
        '''
        at1: str = type1.eqt_bond if type(type1) is AtomType else type1
        at2: str = type2.eqt_bond if type(type2) is AtomType else type2
        bterm = BondTerm(at1, at2, 0)
        try:
            bterm = self.bond_terms[bterm.name]
        except:
            raise FFTermNotFoundError(f'{str(bterm)} not found in FF')

        return bterm

    def get_eqt_for_angle(self, angle):
        '''
        Get the equivalent atom types for angle parameters of the three atoms in a angle.

        Parameters
        ----------
        angle : ~mstools.topology.Angle

        Returns
        -------
        type1 : str
            The equivalent atom type for angle parameters of the first atom in this angle.
        type2 : str
            The equivalent atom type for angle parameters of the second atom in this angle.
        type3 : str
            The equivalent atom type for angle parameters of the third atom in this angle.
        '''
        try:
            at1 = self.atom_types.get(angle.atom1.type).eqt_ang_s
            at2 = self.atom_types.get(angle.atom2.type).eqt_ang_c
            at3 = self.atom_types.get(angle.atom3.type).eqt_ang_s
        except:
            raise FFTermNotFoundError(f'Atom type {angle.atom1.type} or {angle.atom2.type} '
                                      f'or {angle.atom3.type} not found in FF')
        return at1, at2, at3

    def get_angle_term(self, type1, type2, type3):
        '''
        Get the angle term for the angle formed by three atom types.

        If not found, an Exception will be raised.

        Parameters
        ----------
        type1 : AtomType or str
            If type1 is a AtomType object, then will search for its equivalent type for angle side parameters.
            If type1 is a str, then will match exactly this type.
        type2 : AtomType or str
            If type2 is a AtomType object, then will search for its equivalent type for angle center parameters.
            If type2 is a str, them will match exactly this type.
        type3 : AtomType or str
            If type3 is a AtomType object, then will search for its equivalent type for angle side parameters.
            If type3 is a str, them will match exactly this type.

        Returns
        -------
        term : subclass of AngleTerm
        '''
        at1: str = type1.eqt_ang_s if type(type1) is AtomType else type1
        at2: str = type2.eqt_ang_c if type(type2) is AtomType else type2
        at3: str = type3.eqt_ang_s if type(type3) is AtomType else type3
        aterm = AngleTerm(at1, at2, at3, 0)
        try:
            aterm = self.angle_terms[aterm.name]
        except:
            raise FFTermNotFoundError(f'{str(aterm)} not found in FF')

        return aterm

    def get_eqt_for_dihedral(self, dihedral):
        '''
        Get the list of possible equivalent atom types for dihedral parameters of the four atoms in a dihedral.

        Because the wildcard (*) is allowed for the side atoms of a dihedral term,
        a list of possible equivalent atom types is returned.

        Parameters
        ----------
        dihedral : ~mstools.topology.Dihedral

        Returns
        -------
        types : list of tuple of str
            Each element of the list is a tuple of four str (type1, type2, type3, type4),
            which are the equivalent atom types for dihedral parameters of the four atom in this dihedral.
        '''
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

    def get_dihedral_term(self, type1, type2, type3, type4):
        '''
        Get the dihedral term for the dihedral formed by four atom types.

        If not found, an Exception will be raised.

        Parameters
        ----------
        type1 : AtomType or str
            If type1 is a AtomType object, then will search for its equivalent type for dihedral side parameters.
            If type1 is a str, then will match exactly this type.
        type2 : AtomType or str
            If type2 is a AtomType object, then will search for its equivalent type for dihedral center parameters.
            If type2 is a str, them will match exactly this type.
        type3 : AtomType or str
            If type3 is a AtomType object, then will search for its equivalent type for dihedral center parameters.
            If type3 is a str, them will match exactly this type.
        type4 : AtomType or str
            If type4 is a AtomType object, then will search for its equivalent type for dihedral side parameters.
            If type4 is a str, them will match exactly this type.

        Returns
        -------
        term : subclass of DihedralTerm
        '''
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

    def get_eqt_for_improper(self, improper):
        '''
        Get the list of possible equivalent atom types for improper parameters of the four atoms in a improper.

        Because the wildcard (*) is allowed for the side atoms of a improper term,
        a list of possible equivalent atom types is returned.

        Parameters
        ----------
        improper : ~mstools.topology.Improper

        Returns
        -------
        types : list of tuple of str
            Each element of the list is a tuple of four str (type1, type2, type3, type4),
            which are the equivalent atom types for improper parameters of the four atom in this improper.
        '''
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
        '''
        Get the improper term for the improper formed by four atom types.

        If not found, an Exception will be raised.

        Parameters
        ----------
        type1 : AtomType or str
            If type1 is a AtomType object, then will search for its equivalent type for improper center parameters.
            If type1 is a str, then will match exactly this type.
        type2 : AtomType or str
            If type2 is a AtomType object, then will search for its equivalent type for improper side parameters.
            If type2 is a str, them will match exactly this type.
        type3 : AtomType or str
            If type3 is a AtomType object, then will search for its equivalent type for improper side parameters.
            If type3 is a str, them will match exactly this type.
        type4 : AtomType or str
            If type4 is a AtomType object, then will search for its equivalent type for improper side parameters.
            If type4 is a str, them will match exactly this type.

        Returns
        -------
        term : subclass of ImproperTerm
        '''
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

    def get_charge_increment(self, bond):
        '''
        Get the charge increment between the two atoms in a bond.

        Parameters
        ----------
        bond : ~mstools.topology.Bond

        Returns
        -------
        offset : float
            The charge transferred from the second atom to the first atom in this bond.

        '''
        try:
            at1 = self.atom_types.get(bond.atom1.type).eqt_q_inc
            at2 = self.atom_types.get(bond.atom2.type).eqt_q_inc
        except:
            raise FFTermNotFoundError(
                f'AtomType {bond.atom1.type} or {bond.atom2.type} not found in FF')
        qterm = ChargeIncrementTerm(at1, at2, 0)
        try:
            qterm = self.bci_terms[qterm.name]
        except:
            raise FFTermNotFoundError(f'{str(qterm)} not found in FF')

        direction = 1 if at1 == qterm.type1 else -1
        return qterm.value * direction
