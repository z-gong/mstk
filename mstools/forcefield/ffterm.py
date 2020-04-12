from collections import namedtuple


class FFTerm():
    def __init__(self):
        self.name = 'ffterm'

    def __str__(self):
        return '<%s: %s>' % (self.__class__.__name__, self.name)

    def __repr__(self):
        return str(self) + ' at 0x' + str(hex(id(self))[2:].upper())


class AtomType(FFTerm):
    def __init__(self, name, mass=0.0):
        super().__init__()
        self.name = name
        self.mass = mass
        self._charge = 0.0
        # since most of the force fields use LJ-126
        # it's convenient to store LJ-126 parameters here
        # instead of writing every pair as LJTerm
        self._epsilon = 0.0
        self._sigma = 0.0
        self.eqt_vdw = self
        self.eqt_charge = self
        self.eqt_bond_increment = self
        self.eqt_bond = self
        self.eqt_angle_center = self
        self.eqt_angle_side = self
        self.eqt_dihedral_center = self
        self.eqt_dihedral_side = self
        self.eqt_improper_center = self
        self.eqt_improper_side = self

    @property
    def charge(self):
        if self.eqt_charge == self:
            return self._charge
        else:
            return self.eqt_charge._charge

    @charge.setter
    def charge(self, value):
        self._charge = value
        if self.eqt_charge != self:
            print('warning: charge parameter for %s is equivalent to %s, '
                  'you may want to set the charge for %s or change the eqt_charge to self'
                  % (str(self), str(self.eqt_charge), str(self.eqt_charge)))

    @property
    def epsilon(self):
        if self.eqt_vdw == self:
            return self._epsilon
        else:
            return self.eqt_vdw._epsilon

    @epsilon.setter
    def epsilon(self, value):
        self._epsilon = value
        if self.eqt_vdw != self:
            print('warning: epsilon parameter for %s is equivalent to %s, '
                  'you may want to set the epsilon for %s or change the eqt_lj to self'
                  % (str(self), str(self.eqt_charge), str(self.eqt_charge)))

    @property
    def sigma(self):
        if self.eqt_vdw == self:
            return self._sigma
        else:
            return self.eqt_vdw._sigma

    @sigma.setter
    def sigma(self, value):
        self._sigma = value
        if self.eqt_vdw != self:
            print('warning: sigma parameter for %s is equivalent to %s, '
                  'you may want to set the sigma for %s or change the eqt_lj to self'
                  % (str(self), str(self.eqt_charge), str(self.eqt_charge)))

    def __lt__(self, other):
        return self.name < other.name

    def __gt__(self, other):
        return self.name > other.name


class VdwTerm(FFTerm):
    def __init__(self, type1: str, type2: str):
        super().__init__()
        at1, at2 = sorted([type1, type2])
        self.type1 = at1
        self.type2 = at2
        self.name = '%s,%s' % (self.type1, self.type2)

    def __lt__(self, other):
        return [self.type1, self.type2] < [other.type1, other.type2]

    def __gt__(self, other):
        return [self.type1, self.type2] > [other.type1, other.type2]


class ChargeIncrementTerm(FFTerm):
    def __init__(self, type1: str, type2: str, increment: float):
        super().__init__()
        at1, at2 = sorted([type1, type2])
        self.type1 = at1
        self.type2 = at2
        self.increment = increment if at1 == type1 else -increment
        self.name = '%s,%s' % (self.type1, self.type2)

    def __lt__(self, other):
        return [self.type1, self.type2] < [other.type1, other.type2]

    def __gt__(self, other):
        return [self.type1, self.type2] > [other.type1, other.type2]


class BondTerm(FFTerm):
    def __init__(self, type1: str, type2: str, length: float, fixed=False):
        super().__init__()
        at1, at2 = sorted([type1, type2])
        self.type1 = at1
        self.type2 = at2
        self.length = length
        self.fixed = fixed
        self.name = '%s,%s' % (self.type1, self.type2)

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
        self.name = '%s,%s,%s' % (self.type1, self.type2, self.type3)

    def __lt__(self, other):
        return [self.type1, self.type2, self.type3] < [other.type1, other.type2, other.type3]

    def __gt__(self, other):
        return [self.type1, self.type2, self.type3] > [other.type1, other.type2, other.type3]


class DihedralTerm(FFTerm):
    def __init__(self, type1: str, type2: str, type3: str, type4: str):
        super().__init__()
        at1, at2, at3, at4 = min([(type1, type2, type3, type4), (type4, type3, type2, type1)])
        self.type1 = at1
        self.type2 = at2
        self.type3 = at3
        self.type4 = at4
        self.name = '%s,%s,%s,%s' % (self.type1, self.type2, self.type3, self.type4)

    def __lt__(self, other):
        return [self.type1, self.type2, self.type3, self.type4] \
               < [other.type1, other.type2, other.type3, other.type4]

    def __gt__(self, other):
        return [self.type1, self.type2, self.type3, self.type4] \
               > [other.type1, other.type2, other.type3, other.type4]


class ImproperTerm(FFTerm):
    '''
    Center atom is the first, following the convention of GROMACS
    '''

    def __init__(self, type1: str, type2: str, type3: str, type4: str):
        super().__init__()
        at2, at3, at4 = sorted([type2, type3, type4])
        self.type1 = type1
        self.type2 = at2
        self.type3 = at3
        self.type4 = at4
        self.name = '%s,%s,%s,%s' % (self.type1, self.type2, self.type3, self.type4)

    def __lt__(self, other):
        return [self.type1, self.type2, self.type3, self.type4] \
               < [other.type1, other.type2, other.type3, other.type4]

    def __gt__(self, other):
        return [self.type1, self.type2, self.type3, self.type4] \
               > [other.type1, other.type2, other.type3, other.type4]


class LJ126Term(VdwTerm):
    '''
    LJ126 is commonly used, so don't generalize it with LJTerm
    '''
    def __init__(self, type1, type2, epsilon, sigma):
        super().__init__(type1, type2)
        self.epsilon = epsilon
        self.sigma = sigma


class LJTerm(VdwTerm):
    '''
    Mainly used for SDK and SAFT-gamma coarse-grained force field
    '''
    def __init__(self, type1, type2, epsilon, sigma, repulsion, attraction):
        super().__init__(type1, type2)
        self.epsilon = epsilon
        self.sigma = sigma
        self.repulsion = repulsion
        self.attraction = attraction


class HarmonicBondTerm(BondTerm):
    '''
    U = k * (b-b0)^2
    '''

    def __init__(self, type1, type2, length, k, fixed=False):
        super().__init__(type1, type2, length, fixed=fixed)
        self.k = k


class HarmonicAngleTerm(AngleTerm):
    '''
    U = k * (theta-theta0)^2
    '''

    def __init__(self, type1, type2, type3, theta, k, fixed=False):
        super().__init__(type1, type2, type3, theta, fixed=fixed)
        self.k = k


class PeriodicDihedralTerm(DihedralTerm):
    '''
    U = k * (1+cos(n*phi-phi0))
    '''
    Parameter = namedtuple('Parameter', ['multiplicity', 'phi', 'k'])

    def __init__(self, type1, type2, type3, type4):
        super().__init__(type1, type2, type3, type4)
        self.parameters = []

    def add_parameter(self, multiplicity, phi, k):
        if any([p.multiplicity == multiplicity for p in self.parameters]):
            raise Exception('Parameter already exist with multiplicity %i' % multiplicity)
        para = PeriodicDihedralTerm.Parameter(multiplicity=multiplicity, phi=phi, k=k)
        self.parameters.append(para)


class FourierDihedralTerm(DihedralTerm):
    '''
    U = (k1*(1+cos(phi)) + k2*(1-cos(2*phi)) + k3*(1+cos(3*phi)) + k4*(1-cos(4*phi)))
    '''

    def __init__(self, type1, type2, type3, type4, k1, k2, k3, k4):
        super().__init__(type1, type2, type3, type4)
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4


class PeriodicImproperTerm(ImproperTerm):
    '''
    U = k * (1-cos(2*phi))
    '''

    def __init__(self, type1, type2, type3, type4, k):
        super().__init__(type1, type2, type3, type4)
        self.k = k


class HarmonicImproperTerm(ImproperTerm):
    '''
    U = k * phi^2
    '''

    def __init__(self, type1, type2, type3, type4, k):
        super().__init__(type1, type2, type3, type4)
        self.k = k
