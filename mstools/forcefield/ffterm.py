from collections import namedtuple


class FFTerm():
    def __init__(self):
        pass

    def __repr__(self):
        return '<%s: %s>' % (self.__class__.__name__, self.key)

    @property
    def key(self):
        return 'ffterm'


class AtomType(FFTerm):
    def __init__(self, name, mass=0.0):
        super().__init__()
        self.name = name
        self.mass = mass
        self._charge = 0.0
        self._epsilon = 0.0
        self._sigma = 0.0
        self.equivalent_type_lj = self
        self.equivalent_type_charge = self
        self.equivalent_type_bond_increment = self
        self.equivalent_type_bond = self
        self.equivalent_type_angle_center = self
        self.equivalent_type_angle_side = self
        self.equivalent_type_dihedral_center = self
        self.equivalent_type_dihedral_side = self
        self.equivalent_type_improper_center = self
        self.equivalent_type_improper_side = self

    @property
    def charge(self):
        if self.equivalent_type_charge == self:
            return self._charge
        else:
            return self.equivalent_type_charge._charge

    @charge.setter
    def charge(self, value):
        self._charge = value
        if self.equivalent_type_charge != self:
            print('warning: charge parameter for %s is equivalent to %s, '
                  'you may want to set the charge for %s or change the equivalent_type_charge to self'
                  % (self, self.equivalent_type_charge, self.equivalent_type_charge))

    @property
    def epsilon(self):
        if self.equivalent_type_lj == self:
            return self._epsilon
        else:
            return self.equivalent_type_lj._epsilon

    @epsilon.setter
    def epsilon(self, value):
        self._epsilon = value
        if self.equivalent_type_lj != self:
            print('warning: epsilon parameter for %s is equivalent to %s, '
                  'you may want to set the epsilon for %s or change the equivalent_type_lj to self'
                  % (self, self.equivalent_type_lj, self.equivalent_type_lj))

    @property
    def sigma(self):
        if self.equivalent_type_lj == self:
            return self._sigma
        else:
            return self.equivalent_type_lj._sigma

    @sigma.setter
    def sigma(self, value):
        self._sigma = value
        if self.equivalent_type_lj != self:
            print('warning: sigma parameter for %s is equivalent to %s, '
                  'you may want to set the sigma for %s or change the equivalent_type_lj to self'
                  % (self, self.equivalent_type_lj, self.equivalent_type_lj))

    @property
    def key(self):
        return self.name


class BondTerm(FFTerm):
    def __init__(self, atom1, atom2):
        super().__init__()
        at1, at2 = sorted([atom1, atom2])
        self.atom1 = at1
        self.atom2 = at2

    @property
    def key(self):
        return '%s,%s' % (self.atom1, self.atom2)


class AngleTerm(FFTerm):
    def __init__(self, atom1, atom2, atom3):
        super().__init__()
        at1, at3 = sorted([atom1, atom3])
        self.atom1 = at1
        self.atom2 = atom2
        self.atom3 = at3

    @property
    def key(self):
        return '%s,%s,%s' % (self.atom1, self.atom2, self.atom3)


class DihedralTerm(FFTerm):
    def __init__(self, atom1, atom2, atom3, atom4):
        super().__init__()
        at1, at2, at3, at4 = min([(atom1, atom2, atom3, atom4), (atom4, atom3, atom2, atom1)])
        self.atom1 = at1
        self.atom2 = at2
        self.atom3 = at3
        self.atom4 = at4

    @property
    def key(self):
        return '%s,%s,%s,%s' % (self.atom1, self.atom2, self.atom3, self.atom4)


class ImproperTerm(FFTerm):
    '''
    center atom is the first, following the convention of GROMACS
    '''
    def __init__(self, atom1, atom2, atom3, atom4):
        super().__init__()
        at2, at3, at4 = sorted([atom2, atom3, atom4])
        self.atom1 = atom1
        self.atom2 = at2
        self.atom3 = at3
        self.atom4 = at4

    @property
    def key(self):
        return '%s,%s,%s,%s' % (self.atom1, self.atom2, self.atom3, self.atom4)


class NonbondedTerm(FFTerm):
    def __init__(self, atom1, atom2):
        super().__init__()
        at1, at2 = sorted([atom1, atom2])
        self.atom1 = at1
        self.atom2 = at2

    @property
    def key(self):
        return '%s,%s' % (self.atom1, self.atom2)


class ChargeIncrementTerm(FFTerm):
    def __init__(self, atom1, atom2, increment):
        super().__init__()
        at1, at2 = sorted([atom1, atom2])
        self.atom1 = at1
        self.atom2 = at2
        self.increment = increment if at1 == atom1 else -increment

    @property
    def key(self):
        return '%s,%s' % (self.atom1, self.atom2)


class HarmonicBondTerm(BondTerm):
    '''
    U = 0.5 * k * (b-b0)^2
    '''

    def __init__(self, atom1, atom2, length, k):
        super().__init__(atom1, atom2)
        self.length = length
        self.k = k


class HarmonicAngleTerm(AngleTerm):
    '''
    U = 0.5 * k * (theta-theta0)^2
    '''

    def __init__(self, atom1, atom2, atom3, theta, k):
        super().__init__(atom1, atom2, atom3)
        self.theta = theta
        self.k = k


class PeriodicDihedralTerm(DihedralTerm):
    '''
    U = k*(1+cos(n*phi-phi0))
    '''
    Parameter = namedtuple('Parameter', ['multiplicity', 'phi', 'k'])

    def __init__(self, atom1, atom2, atom3, atom4):
        super().__init__(atom1, atom2, atom3, atom4)
        self.parameters = []

    def add_parameter(self, multiplicity, phi, k):
        if any([p.multiplicity == multiplicity for p in self.parameters]):
            raise Exception('Parameter already exist with multiplicity %i' % multiplicity)
        para = PeriodicDihedralTerm.Parameter(multiplicity=multiplicity, phi=phi, k=k)
        self.parameters.append(para)


class FourierDihedralTerm(DihedralTerm):
    '''
    U = 0.5 * (k1*(1+cos(phi)) + k2*(1-cos(2*phi)) + k3*(1+cos(3*phi)) + k4*(1-cos(4*phi)))
    '''

    def __init__(self, atom1, atom2, atom3, atom4, k1, k2, k3, k4):
        super().__init__(atom1, atom2, atom3, atom4)
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4


class PeriodicImproperTerm(ImproperTerm):
    '''
    U = k * (1-cos(2*phi))
    '''

    def __init__(self, atom1, atom2, atom3, atom4, k):
        super().__init__(atom1, atom2, atom3, atom4)
        self.k = k


class HarmonicImproperTerm(ImproperTerm):
    '''
    U = 0.5 * k * phi^2
    '''

    def __init__(self, atom1, atom2, atom3, atom4, k):
        super().__init__(atom1, atom2, atom3, atom4)
        self.k = k


class LJTerm(NonbondedTerm):
    def __init__(self, atom1, atom2, epsilon, sigma, repulsion, attraction):
        super().__init__(atom1, atom2)
        self.epsilon = epsilon
        self.sigma = sigma
        self.repulsion = repulsion
        self.attraction = attraction
