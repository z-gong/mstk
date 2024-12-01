import math
from distutils.util import strtobool
from mstk import logger
from mstk.chem.constant import *
from .errors import *


class FFTermFactory:
    '''
    Factory class for force field terms.

    Mainly used for loading a FFTerm from a line in ZFF file.
    '''

    _class_map = {}

    @staticmethod
    def register(klass, extra_names=None):
        '''
        Register a force field term so that it can be read from or write into a ZFF file.

        Parameters
        ----------
        klass : subclass of FFTerm
        extra_names : list of str
        '''
        FFTermFactory._class_map[klass.__name__] = klass
        FFTermFactory._class_map[klass.get_alias()] = klass  # for short alias in ZFF
        extra_names = extra_names or []
        for name in extra_names:
            FFTermFactory._class_map[name] = klass

    @staticmethod
    def create_from_zfp(name, d):
        '''
        Reconstruct a FFTerm using a dict of attributes read from ZFP xml file.

        Essentially it finds the class for the term, and calls the constructor of the class.
        Every class to be reconstructed should have a class attribute `_term_attrs`.
        This dict records which attributes should be saved into ZFP xml file,
        and it is also used as the arguments for initializing a FFTerm.
        It also tells the data type of these attributes.

        Some force field terms can not be fully built from the constructor.
        Then :func:`FFTerm._from_zfp_extra` should be implemented for those terms,
        and it will be called to do some extra works.

        Parameters
        ----------
        name : str
            Class name of the term to be constructed
        d : dict, [str, str]
            Attributes loaded from ZFP xml element

        Returns
        -------
        term : subclass of FFTerm

        '''
        # must use kwargs instead of args to init the cls,
        # in case there is mistake in the sequence of init arguments
        try:
            cls = FFTermFactory._class_map[name]
        except:
            raise Exception('Invalid force field term: %s' % name)
        kwargs = {}
        adjustables = []
        for attr, func in cls._term_attrs.items():
            str_val = d[attr]
            if str_val.endswith('?'):
                str_val = str_val[:-1]
                if func is float:
                    adjustables.append(attr)
            val = func(str_val)
            if attr in cls._attr_convert:
                val = cls._attr_convert[attr][0](val)
            kwargs[attr] = val
        try:
            term = cls(**kwargs)
        except:
            raise Exception('Cannot init %s with kwargs %s' % (name, str(d)))
        term.adjustables = adjustables
        term._from_zfp_extra(d)
        return term

    @staticmethod
    def create_from_zff(line):
        '''
        Reconstruct a FFTerm using a line read from ZFF text file.

        Essentially it finds the class for the term, and calls the constructor of the class.
        Every class to be reconstructed should have a class attribute `_term_attrs`.
        This dict records which attributes should be saved into ZFF file,
        and it is also used as the arguments for initializing a FFTerm.
        It also tells the data type of these attributes.

        Some force field terms can not be fully built from the constructor.
        Then :func:`FFTerm._from_zff_extra` should be implemented for those terms,
        and it will be called to do some extra works.

        Parameters
        ----------
        line : str
            A line read from ZFF file

        Returns
        -------
        term : subclass of FFTerm

        '''
        # must use kwargs instead of args to init the cls,
        # in case there is mistake in the sequence of init arguments
        words = line.strip().split()
        name = words[0]
        try:
            cls = FFTermFactory._class_map[name]
        except KeyError:
            raise Exception('Invalid force field term: %s' % name)
        kwargs = {}
        adjustables = []
        for (attr, func), str_val in zip(cls._term_attrs.items(), words[1:]):
            if str_val.endswith('?'):
                str_val = str_val[:-1]
                if func is float:
                    adjustables.append(attr)
            val = func(str_val)
            if attr in cls._attr_convert:
                val = cls._attr_convert[attr][0](val)
            kwargs[attr] = val
        try:
            term = cls(**kwargs)
        except:
            raise Exception('Cannot init %s with line: %s' % (name, line))
        term.adjustables = adjustables
        term._from_zff_extra(line)

        return term


class FFTerm:
    '''
    Base class for all force field terms like :class:`AtomType`, :class:`BondTerm`, etc...

    Attributes
    ----------
    version : str
        Version is used mainly in :class:`PPF` format by DFF.
    comments : list of str
        Comments is useful for debugging when generating input files for simulation engines.
    '''

    _term_attrs = {}
    '''
    The attributes that will be saved into or loaded from ZFF file.
    The key is the attribute name. The value is a function to parse the string to object.
    This dict must be overridden by every subclasses that can be saved into or reconstructed from ZFF file.
    '''

    _attr_convert = {}
    '''
    The parameters that will be converted when read/write parameters from/to ZFF file.
    The key is the attribute name. The value is a tuple of two functions.
    The first one convert data from external to internal. The second one do the opposite.
    For example, the angle is in unit of radian internally. But it is in unit of degree in ZFF file.
    Therefore, for AngleTerm, _attr_convert have two functions for attribute 'theta'.
    The first one convert degree to radian, and the second one convert radian to degree.
    This dict must be overridden by subclasses that have different value internally and in ZFF file.
    '''

    @classmethod
    def get_alias(cls):
        '''
        Return a short alias for the name of this class

        It is used to keep the line more compact in ZFF file
        '''
        return cls.__name__.replace('Term', '')

    def __init__(self):
        self.version = None
        self.comments: [str] = []
        self.adjustables: [str] = []  # parameters that are subject to optimization

    def __repr__(self):
        return '<%s: %s>' % (self.__class__.__name__, self.name)

    @property
    def name(self):
        '''
        The name of this force field term.

        Returns
        -------
        name : str
        '''
        raise NotImplementedError('This property should be implemented by subclasses')

    def _from_zfp_extra(self, d):
        '''
        Some force field terms can not be fully built from the constructor.
        This function is called by :func:`FFTermFactory.create_from_zfp` to do some extra works.

        e.g. PeriodicDihedralTerm need to rebuild its own storage of multiple dihedral parameters

        Parameters
        ----------
        d : dict, [str, str]
        '''
        pass

    def to_zff(self):
        '''
        Pack the attributes of a term into a string so that can be saved into a line in ZFF file.

        Every class to be packed should have a class property `_term_attrs`.
        This dict records which attributes should be saved into ZFF file,
        and it is also used as the arguments for initializing a FFTerm.
        It also tells the data type of these attributes.

        Returns
        -------
        line : string
        '''
        line = '%-16s' % self.__class__.get_alias()
        for attr, func in self._term_attrs.items():
            val = getattr(self, attr)
            if attr in self._attr_convert:
                val = self._attr_convert[attr][1](val)
            if func is float:
                if abs(val) >= 1E6:
                    string = ' %8i' % round(val)
                elif abs(val) >= 1E3:
                    string = ' %8.1f' % val
                elif 0 < abs(val) < 1E-3:
                    string = ' %8.1e' % val
                else:
                    string = ' %8.4f' % val
                string += '?' if attr in self.adjustables else ' '
            elif func is int:
                string = ' %2i' % val
            else:
                string = ' %-9s' % val
            line += string

        return line.rstrip()

    def to_zff_header(self):
        '''
        Header string to explain a line in ZFF format

        Returns
        -------
        header : str
        '''
        line = '#%-15s' % self.__class__.get_alias()
        for attr, func in self._term_attrs.items():
            if func is float:
                string = ' %8s ' % attr
            elif func is int:
                string = ' %2i' % attr
            else:
                string = ' %-9s' % attr
            line += string
        return line.rstrip()

    def _from_zff_extra(self, line):
        '''
        Some force field terms can not be fully built from the constructor.
        This function is called by :func:`FFTermFactory.create_from_zff` to do some extra works.

        e.g. PeriodicDihedralTerm need to rebuild its own storage of multiple dihedral parameters

        Parameters
        ----------
        line : str
        '''
        pass

    def evaluate_energy(self, val):
        '''
        Evaluate the energy for a force field term like HarmonicBondTerm, PeriodicDihedralTerm, etc...

        It is mainly for debugging. It is not (and can not be) implemented for all terms.

        Parameters
        ----------
        val : float
            The value of distance, angle, dihedral or improper.

        Returns
        -------
        energy : float
        '''
        raise NotImplementedError('This method is invalid for this term')


class AtomType(FFTerm):
    '''
    AtomType is the most fundamental element in a force field.
    It determines which :class:`VdwTerm`, :class:`BondTerm`, etc...
    will be used to describe the energy functions between a specific set of atoms.

    Equivalent table (EQT) is extensively used in `mstk`.
    It helps decreasing the number of force field terms significantly without losing generality and accuracy,
    thus making the force field extendable.

    Parameters
    ----------
    name : str
        The name of this atom type.
    mass : float , optional
        The mass of this atom type.
    charge : float , optional
        The charge of this atom type.
    eqt_bci : str , optional
        The equivalent type for bond charge increment parameters
    eqt_vdw : str , optional
        The equivalent type for vdW parameters
    eqt_bond : str , optional
        The equivalent type for bond parameters
    eqt_ang_c : str , optional
        The equivalent type for angle parameters if this atom type is the center of an angle
    eqt_ang_s: str , optional
        The equivalent type for angle parameters if this atom type is the side of an angle
    eqt_dih_c : str , optional
        The equivalent type for dihedral parameters if this atom type is the center of a dihedral
    eqt_dih_s: str , optional
        The equivalent type for dihedral parameters if this atom type is the side of a dihedral
    eqt_imp_c : str , optional
        The equivalent type for improper parameters if this atom type is the center of an improper
    eqt_imp_s : str , optional
        The equivalent type for improper parameters if this atom type is the side of an improper
    eqt_polar : str , optional
        The equivalent type for polarization parameters

    Attributes
    ----------
    mass : float
        The mass of this atom type. -1 means not provided in the force field.
        Mass is not mandatory, because finally we are going to take the masses from the :class:`~mstk.topology.Topology`.
        However, masses in force field can be used to assign the masses in topology.
    charge : float
        The charge of this atom type.
        Similar to mass, the charges we are finally going to use are from the :class:`~mstk.topology.Topology`.
        But charges in force field can be used to assign the charges in topology.
    eqt_bci : str
        The equivalent type for bond charge increment parameters
    eqt_vdw : str
        The equivalent type for vdW parameters
    eqt_bond : str
        The equivalent type for bond parameters
    eqt_ang_c : str
        The equivalent type for angle parameters if this atom type is the center of an angle
    eqt_ang_s: str
        The equivalent type for angle parameters if this atom type is the side of an angle
    eqt_dih_c : str
        The equivalent type for dihedral parameters if this atom type is the center of a dihedral
    eqt_dih_s: str
        The equivalent type for dihedral parameters if this atom type is the side of a dihedral
    eqt_imp_c : str
        The equivalent type for improper parameters if this atom type is the center of an improper
    eqt_imp_s : str
        The equivalent type for improper parameters if this atom type is the side of an improper
    eqt_polar : str
        The equivalent type for polarization parameters
    '''
    _term_attrs = {
        'name'     : str,
        'mass'     : float,
        'charge'   : float,
        'eqt_bci'  : str,
        'eqt_vdw'  : str,
        'eqt_bond' : str,
        'eqt_ang_c': str,
        'eqt_ang_s': str,
        'eqt_dih_c': str,
        'eqt_dih_s': str,
        'eqt_imp_c': str,
        'eqt_imp_s': str,
        'eqt_polar': str,
    }

    def __init__(self, name, mass=-1, charge=0,
                 eqt_bci=None, eqt_vdw=None, eqt_bond=None,
                 eqt_ang_c=None, eqt_ang_s=None, eqt_dih_c=None, eqt_dih_s=None,
                 eqt_imp_c=None, eqt_imp_s=None, eqt_polar=None):
        super().__init__()
        # cannot set name directly because name is a property for each FFTerm subclass
        # to allow create_from_zff() work correctly for AtomType, the name must be set with @name.setter
        self._name = name
        # mass is only used for assign mass of topology in case mass not exist in topology
        # -1 means mass unknown
        self.mass = mass
        self.charge = charge
        self.eqt_bci = eqt_bci or name
        self.eqt_vdw = eqt_vdw or name
        self.eqt_bond = eqt_bond or name
        self.eqt_ang_c = eqt_ang_c or name
        self.eqt_ang_s = eqt_ang_s or name
        self.eqt_dih_c = eqt_dih_c or name
        self.eqt_dih_s = eqt_dih_s or name
        self.eqt_imp_c = eqt_imp_c or name
        self.eqt_imp_s = eqt_imp_s or name
        self.eqt_polar = eqt_polar or name

    @property
    def name(self):
        '''
        The name of this atom type.

        :setter: Set the name of this atom type

        Returns
        -------
        name : str
        '''
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def eqt_types(self):
        return [self.eqt_bci,
                self.eqt_vdw,
                self.eqt_bond,
                self.eqt_ang_c,
                self.eqt_ang_s,
                self.eqt_dih_c,
                self.eqt_dih_s,
                self.eqt_imp_c,
                self.eqt_imp_s,
                self.eqt_polar]

    def __lt__(self, other):
        return self.name < other.name

    def __gt__(self, other):
        return self.name > other.name


FFTermFactory.register(AtomType)


class ChargeIncrementTerm(FFTerm):
    '''
    ChargeIncrementTerm describes the offset of charge between two atom types when they are bonded.

    Two atom types and the offset from the second to first should be provided to construct this term.
    e.g. The following example means when h_1 and c_4 are bonded,
    0.06 element charges will be transferred from c_4 to h_1.

    >>> term = ChargeIncrementTerm('h_1', 'c_4', 0.06)

    During the initialization, the two atom types will be sorted by their string.
    The value will be modified to represent the offset from the second to the first atom types after sorting.

    >>> print(term.type1, term.type2, term.value)
        c_4 h_1 -0.06

    Parameters
    ----------
    type1 : str
        The string of first atom type
    type2 : str
        The string of second atom type
    value : float
        The offset from type2 to type1

    Attributes
    ----------
    type1 : str
        The string of first atom type
    type2 : str
        The string of second atom type
    value : float
        The offset from type2 to type1
    '''
    _term_attrs = {
        'type1': str,
        'type2': str,
        'value': float,
    }

    def __init__(self, type1, type2, value):
        super().__init__()
        if type1 == type2 and value != 0:
            raise ChargeIncrementNonZeroError('Non-zero charge increment between atoms of same type')

        at1, at2 = sorted([type1, type2])
        self.type1 = at1
        self.type2 = at2
        self.value = value if at1 == type1 else -value

    @property
    def name(self):
        return '%s,%s' % (self.type1, self.type2)

    def __lt__(self, other):
        return [self.type1, self.type2] < [other.type1, other.type2]

    def __gt__(self, other):
        return [self.type1, self.type2] > [other.type1, other.type2]


FFTermFactory.register(ChargeIncrementTerm)


class VdwTerm(FFTerm):
    '''
    Base class for all vdW terms.

    Subclasses should be implemented for each functional form,
    e.g. :class:`LJ126Term`, :class:`MieTerm`.

    During the initialization, the two atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str

    Attributes
    ----------
    type1 : str
    type2 : str
    '''

    def __init__(self, type1, type2):
        super().__init__()
        at1, at2 = sorted([type1, type2])
        self.type1 = at1
        self.type2 = at2

    @property
    def name(self):
        return '%s,%s' % (self.type1, self.type2)

    def __lt__(self, other):
        return [self.type1, self.type2] < [other.type1, other.type2]

    def __gt__(self, other):
        return [self.type1, self.type2] > [other.type1, other.type2]


class BondTerm(FFTerm):
    '''
    Base class for all bond terms.

    Subclasses should be implemented for each functional form,
    e.g. :class:`HarmonicBondTerm`.

    During the initialization, the two atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str
    length : float
    fixed : bool

    Attributes
    ----------
    type1 : str
    type2 : str
    length : float
    fixed : bool
    '''

    def __init__(self, type1, type2, length, fixed=False):
        super().__init__()
        at1, at2 = sorted([type1, type2])
        self.type1 = at1
        self.type2 = at2
        self.length = length
        self.fixed = fixed

    @property
    def name(self):
        return '%s,%s' % (self.type1, self.type2)

    def __lt__(self, other):
        return [self.type1, self.type2] < [other.type1, other.type2]

    def __gt__(self, other):
        return [self.type1, self.type2] > [other.type1, other.type2]


class AngleTerm(FFTerm):
    '''
    Base class for all angle terms.

    Subclasses should be implemented for each functional form,
    e.g. :class:`HarmonicAngleTerm`, :class:`SDKAngleTerm`.

    During the initialization, the two side atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str
    type3 : str
    theta : float
    fixed : bool

    Attributes
    ----------
    type1 : str
    type2 : str
    type3 : str
    theta : float
    fixed : bool
    '''

    # the angle in ZFP or ZFF file is in unit of degree. one digit is enough
    _attr_convert = {
        'theta': (lambda x: x * DEG2RAD, lambda x: round(x * RAD2DEG, 2))
    }

    def __init__(self, type1, type2, type3, theta, fixed=False):
        super().__init__()
        at1, at3 = sorted([type1, type3])
        self.type1 = at1
        self.type2 = type2
        self.type3 = at3
        self.theta = theta
        self.fixed = fixed

    @property
    def name(self):
        return '%s,%s,%s' % (self.type1, self.type2, self.type3)

    def __lt__(self, other):
        return [self.type1, self.type2, self.type3] < [other.type1, other.type2, other.type3]

    def __gt__(self, other):
        return [self.type1, self.type2, self.type3] > [other.type1, other.type2, other.type3]

    @property
    def is_linear(self):
        '''
        Whether or not this angle term is linear.

        Returns
        -------
        is :  bool
        '''
        return False


class DihedralTerm(FFTerm):
    '''
    Base class for all dihedral terms.

    Subclasses should be implemented for each functional form,
    e.g. :class:`PeriodicDihedralTerm`.

    During the initialization, the sequence of atom types to use `i-j-k-l` or `l-k-j-i`
    will be determined by their string order.

    DihedralTerm allows wildcard(*) for side atoms.
    For sorting purpose, the wildcard will always be the last.
    During force field matching, the terms with wildcard will have lower priority.
    The terms with wildcard will be used only if exact match cannot be found.

    Parameters
    ----------
    type1 : str
    type2 : str
    type3 : str
    type4 : str

    Attributes
    ----------
    type1 : str
    type2 : str
    type3 : str
    type4 : str
    '''

    def __init__(self, type1, type2, type3, type4):
        super().__init__()
        if '*' in [type2, type3]:
            raise Exception('Wildcard not allowed for center atoms in DihedralTerm')
        if [type1, type4].count('*') in [0, 2]:
            at1, at2, at3, at4 = min([(type1, type2, type3, type4), (type4, type3, type2, type1)])
        elif type1 == '*':
            at1, at2, at3, at4 = (type4, type3, type2, type1)
        else:
            at1, at2, at3, at4 = (type1, type2, type3, type4)
        self.type1 = at1
        self.type2 = at2
        self.type3 = at3
        self.type4 = at4

    @property
    def name(self):
        return '%s,%s,%s,%s' % (self.type1, self.type2, self.type3, self.type4)

    def __lt__(self, other):
        return [self.type1, self.type2, self.type3, self.type4] \
               < [other.type1, other.type2, other.type3, other.type4]

    def __gt__(self, other):
        return [self.type1, self.type2, self.type3, self.type4] \
               > [other.type1, other.type2, other.type3, other.type4]

    @property
    def is_zero(self):
        '''
        Whether or not this dihedral terms always give zero energy.

        For linear groups like alkyne and nitrile, dihedral terms are presented in topology.
        But the parameters for these terms usually equal to zero.

        Returns
        -------
        is :  bool
        '''
        return False


class ImproperTerm(FFTerm):
    '''
    Base class for all improper terms.

    Subclasses should be implemented for each functional form,
    e.g. :class:`PeriodicImproperTerm`.

    During the initialization, the three side atom types will be sorted by their string.

    Improper term is usually used to keep 3-coordinated structures in the same plane.
    The first atom is the central atom, following the convention of GROMACS and CHARMM.
    However, the definition of the value of improper i-j-k-l is vague.
    CHARMM defines the value of the improper as the angle between plane i-j-k and j-k-l.
    Whereas OPLS defines the value as the angle between j-k-i and k-i-l.

    ImproperTerm allows wildcard(*) for side atoms.
    During force field matching, the terms with wildcard will have lower priority.
    The terms with wildcard will be used only if exact match cannot be found.

    Parameters
    ----------
    type1 : str
    type2 : str
    type3 : str
    type4 : str

    Attributes
    ----------
    type1 : str
    type2 : str
    type3 : str
    type4 : str
    '''

    def __init__(self, type1, type2, type3, type4):
        super().__init__()
        if type1 == '*':
            raise Exception('Wildcard not allowed for center atoms in ImproperTerm')
        non_wildcard = [t for t in (type2, type3, type4) if t != '*']
        at2, at3, at4 = list(sorted(non_wildcard)) + ['*'] * (3 - len(non_wildcard))
        self.type1 = type1
        self.type2 = at2
        self.type3 = at3
        self.type4 = at4

    @property
    def name(self):
        return '%s,%s,%s,%s' % (self.type1, self.type2, self.type3, self.type4)

    def __lt__(self, other):
        return [self.type1, self.type2, self.type3, self.type4] \
               < [other.type1, other.type2, other.type3, other.type4]

    def __gt__(self, other):
        return [self.type1, self.type2, self.type3, self.type4] \
               > [other.type1, other.type2, other.type3, other.type4]


class PolarTerm(FFTerm):
    '''
    Base class for all polarization terms.

    Subclasses should be implemented for each functional form,
    e.g. :class:`DrudePolarTerm`.

    Parameters
    ----------
    type : str

    Attributes
    ----------
    type : str
    '''

    def __init__(self, type):
        super().__init__()
        self.type = type

    @property
    def name(self):
        return self.type

    def __lt__(self, other):
        return self.type < other.type

    def __gt__(self, other):
        return self.type > other.type


class LJ126Term(VdwTerm):
    '''
    vdW term in LJ-12-6 function form.

    The energy function is

    >>> U = 4*epsilon*((sigma/r)^12-(sigma/r)^6)

    LJ126 is commonly used, therefore we don't generalize it with :class:`MieTerm`.

    During the initialization, the two atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str
    epsilon : float
    sigma : float

    Attributes
    ----------
    type1 : str
    type2 : str
    epsilon : float
    sigma : float

    '''

    _term_attrs = {
        'type1'  : str,
        'type2'  : str,
        'epsilon': float,
        'sigma'  : float,
    }

    def __init__(self, type1, type2, epsilon, sigma):
        super().__init__(type1, type2)
        self.epsilon = epsilon
        self.sigma = sigma

    def evaluate_energy(self, r):
        return 4 * self.epsilon * ((self.sigma / r) ** 12 - (self.sigma / r) ** 6)


FFTermFactory.register(LJ126Term)


class MieTerm(VdwTerm):
    '''
    vdW term with generalized LJ function form.

    The energy function is

    >>> U = C * epsilon * ((sigma/r)^n - (sigma/r)^m)
    >>> C = n/(n-m) * (n/m)^(m/(n-m))
    >>> r_min = (n/m)^(1/(n-m)) * sigma

    This term is mainly used for Coarse-Grained force field.

    During the initialization, the two atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str
    epsilon : float
    sigma : float
    repulsion : float
    attraction : float

    Attributes
    ----------
    type1 : str
    type2 : str
    epsilon : float
    sigma : float
    repulsion : float
    attraction : float
    '''

    _term_attrs = {
        'type1'     : str,
        'type2'     : str,
        'epsilon'   : float,
        'sigma'     : float,
        'repulsion' : float,
        'attraction': float,
    }

    def __init__(self, type1, type2, epsilon, sigma, repulsion, attraction):
        super().__init__(type1, type2)
        self.epsilon = epsilon
        self.sigma = sigma
        self.repulsion = repulsion
        self.attraction = attraction

    def evaluate_energy(self, r):
        return self.factor_energy() * self.epsilon \
               * ((self.sigma / r) ** self.repulsion - (self.sigma / r) ** self.attraction)

    def factor_energy(self):
        '''
        Get the energy pre-factor of the Mie function.

        >>> C = n/(n-m) * (n/m)^(m/(n-m))

        Returns
        -------
        factor : float
        '''
        rep, att = self.repulsion, self.attraction
        return rep / (rep - att) * (rep / att) ** (att / (rep - att))

    def factor_r_min(self):
        '''
        Get the converting factor from sigma to the distance of energy minimum.

        >>> f = (n/m)^(1/(n-m))

        Returns
        -------
        factor : float
        '''
        rep, att = self.repulsion, self.attraction
        return (rep / att) ** (1 / (rep - att))

    @property
    def is_sdk(self):
        '''
        SDK CGFF use LJ-9-6 and LJ-12-4

        Returns
        -------
        is : bool
        '''
        if self.repulsion == 9 and self.attraction == 6:
            return True
        if self.repulsion == 12 and self.attraction == 4:
            return True
        return False


FFTermFactory.register(MieTerm)


class MorseVdwTerm(VdwTerm):
    '''
    vdW term with Morse function form.

    The energy function is

    >>> U = depth*((1-exp(-alpha*(r-r0)))^2-1)

    Note its difference to the MorseBondTerm. MorseVdwTerm has zero energy when two atoms are far away,
    while MorseBondTerm has zero energy at equilibrium distance.

    During the initialization, the two atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str
    depth: float
    alpha: float
    r0: float

    Attributes
    ----------
    type1 : str
    type2 : str
    depth: float
    alpha: float
    r0: float
    '''

    _term_attrs = {
        'type1': str,
        'type2': str,
        'depth': float,
        'alpha': float,
        'r0'   : float,
    }

    def __init__(self, type1, type2, depth, alpha, r0):
        super().__init__(type1, type2)
        self.depth = depth
        self.alpha = alpha
        self.r0 = r0

    def evaluate_energy(self, r):
        return self.depth * ((1 - math.exp(-self.alpha * (r - self.r0))) ** 2 - 1)


FFTermFactory.register(MorseVdwTerm)


class HarmonicBondTerm(BondTerm):
    '''
    Bond term in harmonic function form

    The energy function is

    >>> U = k * (b-b0)^2

    During the initialization, the two atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str
    length : float
    k : float
    fixed : bool

    Attributes
    ----------
    type1 : str
    type2 : str
    length : float
    k : float
    fixed : bool
    '''

    _term_attrs = {
        'type1' : str,
        'type2' : str,
        'length': float,
        'k'     : float,
        'fixed' : lambda x: bool(strtobool(x))
    }

    def __init__(self, type1, type2, length, k, fixed=False):
        super().__init__(type1, type2, length, fixed=fixed)
        self.k = k

    def evaluate_energy(self, val):
        return self.k * (val - self.length) ** 2


FFTermFactory.register(HarmonicBondTerm)


class QuarticBondTerm(BondTerm):
    '''
    Bond term in quartic function form

    The energy function is

    >>> U = k2 * (b-b0)^2 + k3 * (b-b0)^3 + k4 * (b-b0)^4

    During the initialization, the two atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str
    length : float
    k2 : float
    k3 : float
    k4 : float
    fixed : bool

    Attributes
    ----------
    type1 : str
    type2 : str
    length : float
    k2 : float
    k3 : float
    k4 : float
    fixed : bool
    '''

    _term_attrs = {
        'type1' : str,
        'type2' : str,
        'length': float,
        'k2'    : float,
        'k3'    : float,
        'k4'    : float,
        'fixed' : lambda x: bool(strtobool(x))
    }

    def __init__(self, type1, type2, length, k2, k3, k4, fixed=False):
        super().__init__(type1, type2, length, fixed=fixed)
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4

    def evaluate_energy(self, val):
        return self.k2 * (val - self.length) ** 2 + self.k3 * (val - self.length) ** 3 \
               + self.k4 * (val - self.length) ** 4


FFTermFactory.register(QuarticBondTerm)


class MorseBondTerm(BondTerm):
    '''
    Bond term in Morse function form

    The energy function is

    >>> U = depth*(1-exp(-alpha*(r-b0)))^2
    >>> alpha = (k / depth) ** 0.5

    The `k` corresponds to the force constant `k` in HarmonicBondTerm.

    Note its difference to the MorseVdwTerm. MorseVdwTerm has zero energy when two atoms are far away,
    while MorseBondTerm has zero energy at equilibrium distance.

    During the initialization, the two atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str
    length : float
    k : float
    depth : float

    Attributes
    ----------
    type1 : str
    type2 : str
    length : float
    k : float
    depth : float
    '''

    _term_attrs = {
        'type1' : str,
        'type2' : str,
        'length': float,
        'k'     : float,
        'depth' : float,
    }

    def __init__(self, type1, type2, length, k, depth):
        super().__init__(type1, type2, length)
        self.k = k
        self.depth = depth

    def evaluate_energy(self, val):
        alpha = (self.k / self.depth) ** 0.5
        return self.depth * (1 - math.exp(-alpha * (val - self.length))) ** 2


FFTermFactory.register(MorseBondTerm)


class HarmonicAngleTerm(AngleTerm):
    '''
    Angle term in harmonic function form

    The energy function is

    >>> U = k * (theta-theta0)^2

    The angle is in unit of radian, and the force constant is in unit of kJ/mol/rad^2.

    During the initialization, the two side atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str
    type3 : str
    theta : float
    k : float
    fixed : bool

    Attributes
    ----------
    type1 : str
    type2 : str
    type3 : str
    theta : float
    k : float
    fixed : bool
    '''

    _term_attrs = {
        'type1': str,
        'type2': str,
        'type3': str,
        'theta': float,
        'k'    : float,
        'fixed': lambda x: bool(strtobool(x))
    }

    def __init__(self, type1, type2, type3, theta, k, fixed=False):
        super().__init__(type1, type2, type3, theta, fixed=fixed)
        self.k = k

    def evaluate_energy(self, val):
        return self.k * (val - self.theta) ** 2


FFTermFactory.register(HarmonicAngleTerm)


class QuarticAngleTerm(AngleTerm):
    '''
    Angle term in quartic function form

    The energy function is

    >>> U = k2 * (theta-theta0)^2 + k3 * (theta-theta0)^3 + k4 * (theta-theta0)^4

    The angle is in unit of radian, and the force constant is in unit of kJ/mol/rad^2.

    During the initialization, the two side atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str
    type3 : str
    theta : float
    k2 : float
    k3 : float
    k4 : float
    fixed : bool

    Attributes
    ----------
    type1 : str
    type2 : str
    type3 : str
    theta : float
    k2 : float
    k3 : float
    k4 : float
    fixed : bool
    '''

    _term_attrs = {
        'type1': str,
        'type2': str,
        'type3': str,
        'theta': float,
        'k2'   : float,
        'k3'   : float,
        'k4'   : float,
        'fixed': lambda x: bool(strtobool(x))
    }

    def __init__(self, type1, type2, type3, theta, k2, k3, k4, fixed=False):
        super().__init__(type1, type2, type3, theta, fixed=fixed)
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4

    def evaluate_energy(self, val):
        return self.k2 * (val - self.theta) ** 2 + self.k3 * (val - self.theta) ** 3 \
               + self.k4 * (val - self.theta) ** 4


FFTermFactory.register(QuarticAngleTerm)


class HarmonicCosineAngleTerm(AngleTerm):
    '''
    Angle term in harmonic cosine function form

    The energy function is

    >>> U = k / sin(theta0)^2 * (cos(theta) - cos(theta0))^2

    The angle is in unit of radian, and the force constant is in unit of kJ/mol/rad^2.

    During the initialization, the two side atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str
    type3 : str
    theta : float
    k : float
    fixed : bool

    Attributes
    ----------
    type1 : str
    type2 : str
    type3 : str
    theta : float
    k : float
    fixed : bool
    '''

    _term_attrs = {
        'type1': str,
        'type2': str,
        'type3': str,
        'theta': float,
        'k'    : float,
        'fixed': lambda x: bool(strtobool(x))
    }

    def __init__(self, type1, type2, type3, theta, k, fixed=False):
        super().__init__(type1, type2, type3, theta, fixed=fixed)
        self.k = k

    def evaluate_energy(self, val):
        return self.k / math.sin(self.theta) ** 2 * (math.cos(val) - math.cos(self.theta)) ** 2


FFTermFactory.register(HarmonicCosineAngleTerm)


class SDKAngleTerm(AngleTerm):
    '''
    Angle term in SDK function form

    The energy function is

    >>> U = k * (theta-theta0)^2 + LJ96
    >>> LJ96 = 6.75 * epsilon * ((sigma/r)^9 - (sigma/r)^6) + epsilon, r < 1.144714 * sigma

    The angle is in unit of radian, and the force constant is in unit of kJ/mol/rad^2.

    This term itself does not contain parameters `epsilon` and `sigma`.
    Instead, it expects a :class:`MieTerm` between `type1` and `type3` existing in the :class:`ForceField` object.
    If such a MieTerm does not exist, or the repulsion and attraction of this MieTerm do not equal to 9 and 6,
    then the ForceField is invalid and cannot be used to construct a :class:`~mstk.simsys.System`.

    During the initialization, the two side atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str
    type3 : str
    theta : float
    k : float

    Attributes
    ----------
    type1 : str
    type2 : str
    type3 : str
    theta : float
    k : float
    fixed : bool
    '''

    _term_attrs = {
        'type1': str,
        'type2': str,
        'type3': str,
        'theta': float,
        'k'    : float,
    }

    def __init__(self, type1, type2, type3, theta, k):
        super().__init__(type1, type2, type3, theta, fixed=False)
        self.k = k


FFTermFactory.register(SDKAngleTerm)


class LinearAngleTerm(AngleTerm):
    '''
    Linear angle described by the distance between the center atom and its equilibrated position

    For angle i-j-k, the energy function is

    >>> U = k_linear * (r_j - r_j_0)^2
    >>> r_j_0 = a * r_i + (1-a) * r_k
    >>> a = b_jk / (b_ij + b_jk)

    The k_linear can be converted from the k_theta in HarmonicAngleTerm

    >>> k_linear = k_theta * 2 * (b_ij + b_jk)^2 / (b_ij * b_jk)^2

    To be consistent with other angle terms, the k_theta is defined as the force constant, in unit of kJ/mol/rad^2

    This term itself does not contain bond length parameters for i-j and j-k bonds
    Instead, it expects two :class:`BondTerm` between each side atom type and `type2` existing in the :class:`ForceField` object.
    If such a `BondTerm` does not exist, then the ForceField is invalid and cannot be used to construct a :class:`~mstk.simsys.System`.

    During the initialization, the two side atom types will be sorted by their string.

    Parameters
    ----------
    type1 : str
    type2 : str
    type3 : str
    k : float

    Attributes
    ----------
    type1 : str
    type2 : str
    type3 : str
    theta : float
    k : float
    fixed : bool
    '''

    _term_attrs = {
        'type1': str,
        'type2': str,
        'type3': str,
        'k'    : float,
    }

    def __init__(self, type1, type2, type3, k):
        super().__init__(type1, type2, type3, PI, fixed=False)
        self.k = k

    def calc_a_k_linear(self, b12, b23):
        a = b23 / (b12 + b23)
        k_linear = self.k * 2 * (b12 + b23) ** 2 / (b12 * b23) ** 2
        return a, k_linear


FFTermFactory.register(LinearAngleTerm)


class OplsDihedralTerm(DihedralTerm):
    '''
    Dihedral term in OPLS form

    The energy function is

    >>> U = k_1*(1+cos(phi)) + k_2*(1-cos(2*phi)) + k_3*(1+cos(3*phi)) + k_4*(1-cos(4*phi))

    OPLS form is commonly used, therefore we don't generalize it with :class:`PeriodicDihedralTerm`.

    During the initialization, the sequence of atom types to use `i-j-k-l` or `l-k-j-i`
    will be determined by their string order.

    DihedralTerm allows wildcard(*) for side atoms.
    For sorting purpose, the wildcard will always be the last.
    During force field matching, the terms with wildcard will have lower priority.
    The terms with wildcard will be used only if exact match cannot be found.

    Parameters
    ----------
    type1 : str
    type2 : str
    type3 : str
    type4 : str
    k1 : float
    k2 : float
    k3 : float
    k4 : float

    Attributes
    ----------
    type1 : str
    type2 : str
    type3 : str
    type4 : str
    k1 : float
    k2 : float
    k3 : float
    k4 : float
    '''

    _term_attrs = {
        'type1': str,
        'type2': str,
        'type3': str,
        'type4': str,
        'k1'   : float,
        'k2'   : float,
        'k3'   : float,
        'k4'   : float,
    }

    def __init__(self, type1, type2, type3, type4, k1, k2, k3, k4):
        super().__init__(type1, type2, type3, type4)
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4

    def evaluate_energy(self, val):
        energy = 0
        energy += self.k1 * (1 + math.cos(val)) \
                  + self.k2 * (1 - math.cos(2 * val)) \
                  + self.k3 * (1 + math.cos(3 * val)) \
                  + self.k4 * (1 - math.cos(4 * val))
        return energy

    @property
    def is_zero(self):
        return self.k1 == 0 and self.k2 == 0 and self.k3 == 0 and self.k4 == 0

    def to_periodic_term(self):
        '''
        Get the equivalent PeriodicDihedralTerm

        Returns
        -------
        term : PeriodicDihedralTerm
        '''
        term = PeriodicDihedralTerm(self.type1, self.type2, self.type3, self.type4)
        if self.k1 != 0:
            term.add_phase(0., self.k1, 1)
        if self.k2 != 0:
            term.add_phase(PI, self.k2, 2)
        if self.k3 != 0:
            term.add_phase(0., self.k3, 3)
        if self.k4 != 0:
            term.add_phase(PI, self.k4, 4)

        return term


FFTermFactory.register(OplsDihedralTerm)


class PeriodicDihedralTerm(DihedralTerm):
    '''
    Dihedral term in periodic cosine function form

    The energy function is

    >>> U = \sum^n k_n*(1+cos(n*phi-phi0_n))

    During the initialization, the sequence of atom types to use `i-j-k-l` or `l-k-j-i`
    will be determined by their string order.

    DihedralTerm allows wildcard(*) for side atoms.
    For sorting purpose, the wildcard will always be the last.
    During force field matching, the terms with wildcard will have lower priority.
    The terms with wildcard will be used only if exact match cannot be found.

    PeriodicTermDihedral is initialized without phases.
    The parameters for phases should be added by calling :func:`add_phase`.

    Examples
    --------

    >>> term = PeriodicDihedralTerm('h_1', 'c_4', 'c_4', 'h_1')
    >>> term.add_phase(0.0, 0.6485, 1)
    >>> term.add_phase(3.14, 1.0678, 2)
    >>> term.add_phase(0.0, 0.6226, 3)
    >>> for phase in term.phases:
    >>>     print(phase.phi, phase.k, phase.n)
        0.0 0.6485 1
        3.14 1.0678 2
        0.0 0.6226 3
    >>> term = term.to_opls_term()
    >>> print((term.k1, term.k2, term.k3, term.k4))
        (0.6485 1.0678 0.6226 0.0)

    Parameters
    ----------
    type1 : str
    type2 : str
    type3 : str
    type4 : str

    Attributes
    ----------
    type1 : str
    type2 : str
    type3 : str
    type4 : str
    phases : list of namedtuple
        Each element is a namedtuple :attr:`Phase` describing the parameters for each phase.
    '''

    _term_attrs = {
        'type1': str,
        'type2': str,
        'type3': str,
        'type4': str,
    }

    class Phase:
        def __init__(self, phi, k, n):
            if phi < -PI or phi > PI:
                logger.error(f'phi should be in the range of -PI to PI')
            self.phi = phi
            self.k = k
            self.n = n

    def __init__(self, type1, type2, type3, type4):
        super().__init__(type1, type2, type3, type4)
        self.phases = []

    def add_phase(self, phi, k, n):
        '''
        Add parameters for a phase to this dihedral term.

        Multiplicity must be a positive integer. Otherwise, an Exception will be raised.
        If there already exists parameter for this multiplicity, an Exception will be raised.

        Parameters
        ----------
        phi : float
            Phase shift for this multiplicity, in unit of radian
        k : float
            Force constant
        n : int
            Multiplicity
        '''
        if not isinstance(n, int) or n < 1:
            raise Exception('Multiplicity should be a positive integer: %s' % self.name)
        for para in self.phases:
            if para.n == n:
                raise Exception('Duplicated multiplicity: %s' % self.name)
        self.phases.append(PeriodicDihedralTerm.Phase(phi=phi, k=k, n=n))
        self.phases.sort(key=lambda x: x.n)

    def _from_zfp_extra(self, d: {str: str}):
        for key in d.keys():
            if key.startswith('phi_'):
                n = int(key.split('_')[-1])
                phi = float(d['phi_%i' % n]) * DEG2RAD  # the angle in ZFP file is in unit of degree
                k = float(d['k_%i' % n])
                if k != 0:
                    self.add_phase(phi, k, n)

    def to_zff(self):
        line = '%-16s %-9s %-9s %-9s %-9s' % (self.__class__.get_alias(),
                                              self.type1, self.type2, self.type3, self.type4)
        for phase in self.phases:
            # the angle in ZFP or ZFF file is in unit of degree
            line += ' %9.1f %9.4f %2i' % (phase.phi * RAD2DEG, phase.k, phase.n)
        return line

    def _from_zff_extra(self, line):
        words = line.strip().split()
        str_vals = words[5:]
        if len(str_vals) % 3 != 0:
            raise Exception(f'Invalid parameters for {self.__class__.__name__}')
        for i in range(len(str_vals) // 3):
            phi = float(str_vals[i * 3]) * DEG2RAD  # the angle in ZFF file is in unit of degree
            k = float(str_vals[i * 3 + 1])
            n = int(str_vals[i * 3 + 2])
            if k != 0:
                self.add_phase(phi, k, n)

    def evaluate_energy(self, val):
        energy = 0
        for para in self.phases:
            energy += para.k * (1 + math.cos(para.n * val - para.phi))
        return energy

    @property
    def is_zero(self):
        for para in self.phases:
            if para.k != 0:
                return False
        return True

    def to_opls_term(self):
        '''
        Get the equivalent OplsDihedralTerm.

        In OPLS convention, the largest multiplicity equals to 4.
        The phi_0 for multiplicity (1, 2, 3, 4) equal to (0, PI, 0, PI), respectively.
        If it does not follow OPLS convention, an Exception will be raised.

        Returns
        -------
        term : OplsDihedralTerm
        '''
        phases = [p for p in self.phases if p.k != 0]
        for phase in phases:
            if phase.n > 4:
                raise Exception(f'{str(self)} does not follow OPLS convention, n > 4')
        for phase in phases:
            if phase.phi != 0 and abs(phase.phi - PI) > 0.01:
                raise Exception(f'{str(self)} does not follow OPLS convention, phi != 0 or PI')

        k1 = k2 = k3 = k4 = 0.0
        for phase in phases:
            if phase.n == 1:
                if phase.phi == 0:
                    k1 += phase.k
                else:
                    k1 -= phase.k
            elif phase.n == 2:
                if abs(phase.phi - PI) < 0.01:
                    k2 += phase.k
                else:
                    k2 -= phase.k
            elif phase.n == 3:
                if phase.phi == 0:
                    k3 += phase.k
                else:
                    k3 -= phase.k
            elif phase.n == 4:
                if abs(phase.phi - PI) < 0.01:
                    k4 += phase.k
                else:
                    k4 -= phase.k

        return OplsDihedralTerm(self.type1, self.type2, self.type3, self.type4, k1, k2, k3, k4)


FFTermFactory.register(PeriodicDihedralTerm)


class OplsImproperTerm(ImproperTerm):
    '''
    Improper term in OPLS cosine function form

    The energy function is

    >>> U = k * (1-cos(2*phi))

    During the initialization, the three side atom types will be sorted by their string.

    ImproperTerm allows wildcard(*) for side atoms.
    During force field matching, the terms with wildcard will have lower priority.
    The terms with wildcard will be used only if exact match cannot be found.

    Care must be taken when exporting OplsImproperTerm.
    For improper a1-a2-a3-a4, a1 is the center atom.
    In OPLS convention, the value of the improper is defined as the angle between plane a2-a3-a1 and a3-a1-a4.

    Parameters
    ----------
    type1 : str
    type2 : str
    type3 : str
    type4 : str
    k : float

    Attributes
    ----------
    type1 : str
    type2 : str
    type3 : str
    type4 : str
    k : float
    '''

    _term_attrs = {
        'type1': str,
        'type2': str,
        'type3': str,
        'type4': str,
        'k'    : float,
    }

    def __init__(self, type1, type2, type3, type4, k):
        super().__init__(type1, type2, type3, type4)
        self.k = k

    def evaluate_energy(self, val):
        return self.k * (1 - math.cos(2 * val * DEG2RAD))


FFTermFactory.register(OplsImproperTerm)


class HarmonicImproperTerm(ImproperTerm):
    '''
    Improper term in harmonic function form. Mainly used by CHARMM.

    The energy function is

    >>> U = k * (phi-phi0)^2

    The improper is in unit of radian, and the force constant is in unit of kJ/mol/rad^2.

    During the initialization, the three side atom types will be sorted by their string.

    ImproperTerm allows wildcard(*) for side atoms.
    During force field matching, the terms with wildcard will have lower priority.
    The terms with wildcard will be used only if exact match cannot be found.

    Care must be taken when exporting HarmonicImproperTerm.
    For improper a1-a2-a3-a4, a1 is the center atom.
    In CHARMM convention, the improper is defined as the angle between plane a1-a2-a3 and a2-a3-a4.

    Parameters
    ----------
    type1 : str
    type2 : str
    type3 : str
    type4 : str
    phi : float
    k : float

    Attributes
    ----------
    type1 : str
    type2 : str
    type3 : str
    type4 : str
    phi : float
    k : float
    '''

    _term_attrs = {
        'type1': str,
        'type2': str,
        'type3': str,
        'type4': str,
        'phi'  : float,
        'k'    : float,
    }

    # the angle in ZFP or ZFF file is in unit of degree. one digit is enough
    _attr_convert = {
        'phi': (lambda x: x * DEG2RAD, lambda x: round(x * RAD2DEG, 2))
    }

    def __init__(self, type1, type2, type3, type4, phi, k):
        super().__init__(type1, type2, type3, type4)
        self.phi = phi
        self.k = k

    def evaluate_energy(self, val):
        return self.k * ((val - self.phi) * DEG2RAD) ** 2


FFTermFactory.register(HarmonicImproperTerm)


class DrudePolarTerm(PolarTerm):
    '''
    Polarization term described by Drude induced dipole.

    Currently, only the isotropic Drude polarization is implemented.
    The energy function for polarization is

    >>> E = k * d^2

    The charge assigned to the Drude particle is

    >>> q = - sqrt(4 * pi * eps_0 * (2k) * alpha)

    If this atom is bonded with hydrogen atoms, the argument `merge_alpha_H` determines the polarizability of neighbour H atoms.
    The alpha of these H atoms will be merged into the this atom.
    Therefore, DrudePolarTerm for H atom types should not be provided explicitly. Otherwise, there will be double count.
    The H atoms are identified by check the :attr:`~mstk.topology.Atom.symbol` attribute.
    An atom is considered to be hydrogen if symbol equals 'H'.

    Parameters
    ----------
    type : str
    alpha : float
    thole : float
    k : float
    mass : float
    merge_alpha_H : float

    Attributes
    ----------
    type : str
    alpha : float
    thole : float
    k : float
    mass : float
    merge_alpha_H : float
    '''

    _term_attrs = {
        'type'         : str,
        'alpha'        : float,
        'thole'        : float,
        'k'            : float,
        'mass'         : float,
        'merge_alpha_H': float,
    }

    def __init__(self, type: str, alpha, thole, k=4184 / 2 * 100, mass=0.4, merge_alpha_H=0.0):
        super().__init__(type)
        self.alpha = alpha  # nm^3
        self.thole = thole
        self.k = k  # kJ/mol/nm^2
        self.mass = mass
        self.merge_alpha_H = merge_alpha_H

    def to_zff(self):
        line = '%-16s %-9s %8.3e  %8.4f  %8.1f  %8.4f  %8.3e' % (self.__class__.get_alias(),
                                                                 self.type, self.alpha, self.thole, self.k, self.mass,
                                                                 self.merge_alpha_H)
        return line

    def get_charge(self, alpha=None):
        '''
        Get the charge offset between parent atom and Drude particle.

        The charge assigned on the Drude particle is the negative of this charge offset.

        Parameters
        ----------
        alpha : float
            The polarizability of parent atom in unit of nm^3.

        Returns
        -------
        charge : float
            The charge offset from Drude particle to parent atom.
        '''
        if alpha is None:
            alpha = self.alpha
        factor = math.sqrt(1E-6 / AVOGADRO) / ELEMENTARY_CHARGE
        return math.sqrt(4 * PI * VACUUM_PERMITTIVITY * (2 * self.k) * alpha) * factor


FFTermFactory.register(DrudePolarTerm)


class VirtualSiteTerm(FFTerm):
    '''
    Base class for all virtual site terms.

    Subclasses should be implemented for each functional form,
    e.g. :class:`TIP4PSiteTerm`.

    Parameters
    ----------
    type : str

    Attributes
    ----------
    type : str
    '''

    def __init__(self, type):
        super().__init__()
        self.type = type

    @property
    def name(self):
        return self.type

    def __lt__(self, other):
        return self.type < other.type

    def __gt__(self, other):
        return self.type > other.type


class TIP4PSiteTerm(VirtualSiteTerm):
    '''
    TIP4P style virtual site term.

    Parameters
    ----------
    type : str
    type_O : str
    type_H : str
    d : float

    Attributes
    ----------
    type : str
    type_O : str
    type_H : str
    d : float
    '''

    _term_attrs = {
        'type'  : str,
        'type_O': str,
        'type_H': str,
        'd'     : float,
    }

    def __init__(self, type: str, type_O: str, type_H: str, d: float):
        super().__init__(type)
        self.type_O = type_O
        self.type_H = type_H
        self.d = d


FFTermFactory.register(TIP4PSiteTerm)
