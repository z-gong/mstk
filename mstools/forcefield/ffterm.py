from collections import namedtuple
import math
from distutils.util import strtobool
from ..constant import *


class FFTerm():
    def __init__(self):
        self.version = None

    @property
    def name(self):
        raise NotImplementedError('This property should be implemented by subclasses')

    def __str__(self):
        return '<%s: %s>' % (self.__class__.__name__, self.name)

    def __repr__(self):
        return str(self) + ' at 0x' + str(hex(id(self))[2:].upper())

    def to_zfp(self) -> {str: str}:
        d = {}
        for attr, func in self.zfp_attrs.items():
            d[attr] = str(getattr(self, attr))
        d.update(self.to_zfp_extra())
        return d

    def to_zfp_extra(self) -> {str: str}:
        return {}

    @classmethod
    def from_zfp(cls, d: {str: str}):
        term = cls.singleton()
        for attr, func in term.zfp_attrs.items():
            setattr(term, attr, func(d[attr]))
        term.from_zfp_extra(d)
        return term

    def from_zfp_extra(self, d: {str: str}):
        pass

    @staticmethod
    def singleton():
        raise NotImplementedError('This method should be implemented by subclasses')

    @property
    def zfp_attrs(self):
        raise NotImplementedError('This property should be implemented by subclasses')


class AtomType(FFTerm):
    def __init__(self, name):
        super().__init__()
        self._name = name
        # mass is only used for assign mass of topology in case mass not exist in topology
        # mass = -1 means unknown
        self.mass = -1
        self.charge = 0.
        self.eqt_vdw = name
        self.eqt_q_inc = name
        self.eqt_bond = name
        self.eqt_ang_c = name
        self.eqt_ang_s = name
        self.eqt_dih_c = name
        self.eqt_dih_s = name
        self.eqt_imp_c = name
        self.eqt_imp_s = name
        self.eqt_polar = name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        # this is a trick allowing from_zfp() to work correctly
        self._name = value

    def __lt__(self, other):
        return self.name < other.name

    def __gt__(self, other):
        return self.name > other.name

    @staticmethod
    def singleton():
        return AtomType('')

    @property
    def zfp_attrs(self):
        return {
            'name'     : str,
            'mass'     : float,
            'charge'   : float,
            'eqt_vdw'  : str,
            'eqt_bond' : str,
            'eqt_q_inc': str,
            'eqt_ang_c': str,
            'eqt_ang_s': str,
            'eqt_dih_c': str,
            'eqt_dih_s': str,
            'eqt_imp_c': str,
            'eqt_imp_s': str,
            'eqt_polar': str,
        }


class ChargeIncrementTerm(FFTerm):
    def __init__(self, type1: str, type2: str, value: float):
        super().__init__()
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

    @staticmethod
    def singleton():
        return ChargeIncrementTerm('', '', 0.0)

    @property
    def zfp_attrs(self):
        return {
            'type1': str,
            'type2': str,
            'value': float,
        }


class VdwTerm(FFTerm):
    def __init__(self, type1: str, type2: str):
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

    def evaluate_energy(self, r):
        raise NotImplementedError('This method should be implemented by subclasses')


class BondTerm(FFTerm):
    def __init__(self, type1: str, type2: str, length: float, fixed=False):
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
    def __init__(self, type1: str, type2: str, type3: str, theta: float, fixed=False):
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


class DihedralTerm(FFTerm):
    '''
    DihedralTerm allows wildcard(*) for side atoms
    '''

    def __init__(self, type1: str, type2: str, type3: str, type4: str):
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


class ImproperTerm(FFTerm):
    '''
    ImproperTerm allows wildcard(*) for side atoms
    Center atom is the first, following the convention of GROMACS
    '''

    def __init__(self, type1: str, type2: str, type3: str, type4: str):
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


class PolarizableTerm(FFTerm):
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
    LJ126 is commonly used, so don't generalize it with MieTerm
    U = 4*epsilon*((sigma/r)^12-(sigma/r)^6)
    '''

    def __init__(self, type1, type2, epsilon, sigma):
        super().__init__(type1, type2)
        self.epsilon = epsilon
        self.sigma = sigma

    @staticmethod
    def singleton():
        return LJ126Term('', '', 0., 0.)

    @property
    def zfp_attrs(self):
        return {
            'type1'  : str,
            'type2'  : str,
            'epsilon': float,
            'sigma'  : float,
        }

    def evaluate_energy(self, r):
        return 4 * self.epsilon * ((self.sigma / r) ** 12 - (self.sigma / r) ** 6)


class MieTerm(VdwTerm):
    '''
    Mainly used for SDK and SAFT-gamma Coarse-Grained force field
    U = C * epsilon * ((sigma/r)^n - (sigma/r)^m)
    C = n/(n-m) * (n/m)^(m/(n-m))
    r_min = (n/m)^(1/(n-m)) * sigma
    '''

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
        rep, att = self.repulsion, self.attraction
        return rep / (rep - att) * (rep / att) ** (att / (rep - att))

    def factor_r_min(self):
        rep, att = self.repulsion, self.attraction
        return (rep / att) ** (1 / (rep - att))

    @staticmethod
    def singleton():
        return MieTerm('', '', 0., 0., 0., 0.)

    @property
    def zfp_attrs(self):
        return {
            'type1'     : str,
            'type2'     : str,
            'epsilon'   : float,
            'sigma'     : float,
            'repulsion' : float,
            'attraction': float,
        }


class HarmonicBondTerm(BondTerm):
    '''
    U = k * (b-b0)^2
    '''

    def __init__(self, type1, type2, length, k, fixed=False):
        super().__init__(type1, type2, length, fixed=fixed)
        self.k = k

    @staticmethod
    def singleton():
        return HarmonicBondTerm('', '', 0., 0.)

    @property
    def zfp_attrs(self):
        return {
            'type1' : str,
            'type2' : str,
            'length': float,
            'k'     : float,
            'fixed' : lambda x: bool(strtobool(x))
        }


class HarmonicAngleTerm(AngleTerm):
    '''
    U = k * (theta-theta0)^2
    '''

    def __init__(self, type1, type2, type3, theta, k, fixed=False):
        super().__init__(type1, type2, type3, theta, fixed=fixed)
        self.k = k

    @staticmethod
    def singleton():
        return HarmonicAngleTerm('', '', '', 0., 0.)

    @property
    def zfp_attrs(self):
        return {
            'type1': str,
            'type2': str,
            'type3': str,
            'theta': float,
            'k'    : float,
            'fixed' : lambda x: bool(strtobool(x))
        }


class SDKAngleTerm(AngleTerm):
    '''
    U = k * (theta-theta0)^2 + LJ96
    LJ96 = 6.75 * epsilon * ((sigma/r)^9 - (sigma/r)^6) + epsilon, r < 1.144714 * sigma
    '''

    def __init__(self, type1, type2, type3, theta, k, epsilon, sigma):
        super().__init__(type1, type2, type3, theta, fixed=False)
        self.k = k
        self.epsilon = epsilon
        self.sigma = sigma

    @staticmethod
    def singleton():
        return SDKAngleTerm('', '', '', 0., 0., 0., 0.)

    @property
    def zfp_attrs(self):
        return {
            'type1'  : str,
            'type2'  : str,
            'type3'  : str,
            'theta'  : float,
            'k'      : float,
            'epsilon': float,
            'sigma'  : float,
        }


class PeriodicDihedralTerm(DihedralTerm):
    '''
    U = k * (1+cos(n*phi-phi0_n))
    '''

    Parameter = namedtuple('Parameter', ('phi', 'k', 'n'))

    def __init__(self, type1, type2, type3, type4, parameters=[]):
        super().__init__(type1, type2, type3, type4)
        self.parameters = []
        for para in parameters:
            if len(para) != 3:
                raise Exception(
                    'Three values should be provided for each multiplicity: %s' % self.name)
            self.add_parameter(*para)

    def add_parameter(self, phi, k, n):
        if not isinstance(n, int) or n < 1:
            raise Exception('Multiplicity should be a positive integer: %s' % self.name)
        for para in self.parameters:
            if para.n == n:
                raise Exception('Duplicated multiplicity: %s' % self.name)
        self.parameters.append(self.Parameter(phi=phi, k=k, n=n))

    @staticmethod
    def singleton():
        return PeriodicDihedralTerm('', '', '', '')

    @property
    def zfp_attrs(self):
        return {
            'type1': str,
            'type2': str,
            'type3': str,
            'type4': str,
        }

    def to_zfp_extra(self) -> {str: str}:
        d = {}
        for para in self.parameters:
            d.update({
                'phi_%i' % para.n: '%.1f' % para.phi,
                'k_%i' % para.n  : '%.4f' % para.k,
            })
        return d

    def from_zfp_extra(self, d: {str: str}):
        for key in d.keys():
            if key.startswith('phi_'):
                n = int(key.split('_')[-1])
                phi = float(d['phi_%i' % n])
                k = float(d['k_%i' % n])
                self.add_parameter(phi, k, n)


class OplsDihedralTerm(DihedralTerm):
    '''
    OplsDihedral is a more convenient version of PeriodicDihedral
    U = (k1*(1+cos(phi)) + k2*(1-cos(2*phi)) + k3*(1+cos(3*phi)) + k4*(1-cos(4*phi)))
    '''

    def __init__(self, type1, type2, type3, type4, k1, k2, k3, k4):
        super().__init__(type1, type2, type3, type4)
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4

    @staticmethod
    def singleton():
        return OplsDihedralTerm('', '', '', '', 0., 0., 0., 0.)

    @property
    def zfp_attrs(self):
        return {
            'type1': str,
            'type2': str,
            'type3': str,
            'type4': str,
            'k1'   : float,
            'k2'   : float,
            'k3'   : float,
            'k4'   : float,
        }


class PeriodicImproperTerm(ImproperTerm):
    '''
    Normally for keeping 3-coordinated structures in the same plane, in which case xi0=180
    U = k * (1+cos(2*phi-phi0))
    '''

    def __init__(self, type1, type2, type3, type4, phi, k):
        super().__init__(type1, type2, type3, type4)
        self.phi = phi
        self.k = k

    @staticmethod
    def singleton():
        return PeriodicImproperTerm('', '', '', '', 0., 0.)

    @property
    def zfp_attrs(self):
        return {
            'type1': str,
            'type2': str,
            'type3': str,
            'type4': str,
            'phi'  : float,
            'k'    : float,
        }


class HarmonicImproperTerm(ImproperTerm):
    '''
    U = k * (phi-phi0)^2
    '''

    def __init__(self, type1, type2, type3, type4, phi, k):
        super().__init__(type1, type2, type3, type4)
        self.phi = phi
        self.k = k

    @staticmethod
    def singleton():
        return HarmonicImproperTerm('', '', '', '', 0., 0.)

    @property
    def zfp_attrs(self):
        return {
            'type1': str,
            'type2': str,
            'type3': str,
            'type4': str,
            'phi'  : float,
            'k'    : float,
        }


class DrudeTerm(PolarizableTerm):
    '''
    E = k * d^2
    q = sqrt(4 * pi * eps_0 * (2k) * alpha)
    TODO Implement anisotropic Drude polarization
    '''

    def __init__(self, type: str, alpha, thole, merge_alpha_H=0.0):
        super().__init__(type)
        self.alpha = alpha  # nm^3
        self.thole = thole
        self.k = 4184 / 2 * 100  # kJ/mol/nm^2
        self.mass = 0.4
        self.merge_alpha_H = merge_alpha_H

    @staticmethod
    def singleton():
        return DrudeTerm('', 0., 0.)

    @property
    def zfp_attrs(self):
        return {
            'type'         : str,
            'alpha'        : float,
            'thole'        : float,
            'k'            : float,
            'mass'         : float,
            'merge_alpha_H': float,
        }

    def get_charge(self, alpha=None):
        if alpha is None:
            alpha = self.alpha
        factor = math.sqrt(1E-6 / AVOGADRO) / ELEMENTARY_CHARGE
        return math.sqrt(4 * PI * VACUUM_PERMITTIVITY * (2 * self.k) * alpha) * factor
