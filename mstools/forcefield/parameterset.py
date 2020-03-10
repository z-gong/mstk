import math
from .ffterm import *


class ParameterSet():
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
        self.lj_mixing_rule = ParameterSet.LJ_MIXING_NONE
        self.scale_14_vdw = 1.0
        self.scale_14_coulomb = 1.0

    @property
    def use_only_lj126_for_vdw(self):
        '''
        Some logic can be simplified if only LJ126 is used for vdW interactions
        '''
        for vdw in self.vdw_terms.values():
            if type(vdw) != LJ126Term:
                return False
        return True

    def generate_pairwise_lj126_terms(self):
        '''
        Generate pairwise i-j LJ126 terms using combination rule
        The ones already exist will not be touched
        '''
        if self.lj_mixing_rule == self.LJ_MIXING_NONE:
            raise Exception('Mixing rule for LJ126 parameters has been set to LJ_MIXING_NONE')

        n_atom_types = len(self.atom_types)
        for i in range(n_atom_types):
            for j in range(i, n_atom_types):
                atype1 = list(self.atom_types.values())[i].equivalent_type_vdw
                atype2 = list(self.atom_types.values())[j].equivalent_type_vdw
                lj = LJ126Term(atype1.name, atype2.name, 0.0, 1.0)
                if lj.name in self.vdw_terms:
                    continue
                if self.lj_mixing_rule == self.LJ_MIXING_LB:
                    lj.epsilon = math.sqrt(atype1.epsilon * atype2.epsilon)
                    lj.sigma = (atype1.sigma + atype2.sigma) / 2
                elif self.lj_mixing_rule == self.LJ_MIXING_GEOMETRIC:
                    lj.epsilon = math.sqrt(atype1.epsilon * atype2.epsilon)
                    lj.sigma = math.sqrt(atype1.sigma * atype2.sigma)
                else:
                    raise Exception('Unknown mixing rule for LJ126 parameters')
                self.vdw_terms[lj.name] = lj

    def get_vdw_term(self, type1, type2):
        # TODO
        return VdwTerm(type1, type2)


class FFToolParameterSet(ParameterSet):
    def __init__(self, file, mode='r'):
        super().__init__()
        self.lj_mixing_rule = self.LJ_MIXING_GEOMETRIC
        self.scale_14_vdw = 0.5
        self.scale_14_coulomb = 0.5

        self._file = open(file, mode)
        if mode == 'r':
            self.parse()
        elif mode == 'w':
            raise Exception('Writing is not supported for FFToolParameterSet')

    def __del__(self):
        self.close()

    def parse(self):
        self._eq_atom_types = {}  # add a temporary attribute to match equivalent type later
        for line in self._file:
            if line.startswith('#') or line.strip() == '':
                continue

            if line.lower().startswith('atom'):
                section = 'atoms'
                continue
            elif line.lower().startswith('bond'):
                section = 'bonds'
                continue
            elif line.lower().startswith('angl'):
                section = 'angles'
                continue
            elif line.lower().startswith('dihe'):
                section = 'dihedrals'
                continue
            elif line.lower().startswith('impro'):
                section = 'impropers'
                continue

            words = line.strip().split()
            if section == 'atoms':
                self.parse_atom(words)
            elif section == 'bonds':
                self.parse_bond(words)
            elif section == 'angles':
                self.parse_angle(words)
            elif section == 'dihedrals':
                self.parse_dihedral(words)
            elif section == 'impropers':
                self.parse_improper(words)

        self.setup_equivalent_types()

    def parse_atom(self, words):
        atype_eq = AtomType(words[1])
        # add a temporary attribute to match equivalent type later
        atype_eq._eq_type_name = atype_eq.name
        if atype_eq.name not in self._eq_atom_types.keys():
            # add it to the temporary dict, and unique ones will be moved to self.atom_types later
            self._eq_atom_types[atype_eq.name] = atype_eq

        atype = AtomType(words[0], mass=float(words[2]))
        # add a temporary attribute to match equivalent type later
        atype._eq_type_name = words[1]
        atype.charge = float(words[3])
        if words[4] == 'lj':
            atype.sigma = float(words[5]) / 10  # convert from A to nm
            atype.epsilon = float(words[6])
        elif words[4] == 'lj96':
            sigma = float(words[5]) / 10  # convert from A to nm
            epsilon = float(words[6])
            lj = LJTerm(atype.name, atype.name, epsilon, sigma, 9, 6)
            self.vdw_terms[lj.name] = lj
        elif words[4] == 'lj124':
            sigma = float(words[5]) / 10  # convert from A to nm
            epsilon = float(words[6])
            lj = LJTerm(atype.name, atype.name, epsilon, sigma, 12, 4)
            self.vdw_terms[lj.name] = lj
        else:
            raise Exception('Unsupported vdw potential form: %s' % (words[4]))

        if atype.name not in self.atom_types.keys():
            self.atom_types[atype.name] = atype
        else:
            raise Exception('Duplicated atom type: %s' % atype)

    def setup_equivalent_types(self):
        for name, atype in self._eq_atom_types.items():
            if name not in self.atom_types:
                self.atom_types[name] = atype

        for atype in self.atom_types.values():
            atype.equivalent_type_bond = self.atom_types[atype._eq_type_name]
            atype.equivalent_type_angle_center = self.atom_types[atype._eq_type_name]
            atype.equivalent_type_angle_side = self.atom_types[atype._eq_type_name]
            atype.equivalent_type_dihedral_center = self.atom_types[atype._eq_type_name]
            atype.equivalent_type_dihedral_side = self.atom_types[atype._eq_type_name]
            atype.equivalent_type_improper_center = self.atom_types[atype._eq_type_name]
            atype.equivalent_type_improper_side = self.atom_types[atype._eq_type_name]

    def parse_bond(self, words):
        if words[2] not in ['harm', 'cons']:
            raise Exception('Unsupported bond function: %s' % (words[2]))
        bond = HarmonicBondTerm(words[0], words[1], length=float(words[3]) / 10, k=float(words[4]),
                                fixed=(words[2] == 'cons'))
        if bond.name in self.bond_terms.keys():
            raise Exception('Duplicated bond term: %s' % self)
        self.bond_terms[bond.name] = bond

    def parse_angle(self, words):
        if words[3] not in ['harm']:
            raise Exception('Unsupported angle function: %s' % (words[3]))
        angle = HarmonicAngleTerm(words[0], words[1], words[2], theta=float(words[4]), k=float(words[5]))
        if angle.name in self.angle_terms.keys():
            raise Exception('Duplicated angle term: %s' % self)
        self.angle_terms[angle.name] = angle

    def parse_dihedral(self, words):
        if words[4] not in ['opls']:
            raise Exception('Unsupported dihedral function: %s' % (words[4]))
        dihedral = FourierDihedralTerm(words[0], words[1], words[2], words[3],
                                       k1=float(words[5]), k2=float(words[6]),
                                       k3=float(words[7]), k4=float(words[8])
                                       )
        if dihedral.name in self.dihedral_terms.keys():
            raise Exception('Duplicated dihedral term: %s' % self)
        self.dihedral_terms[dihedral.name] = dihedral

    def parse_improper(self, words):
        '''
        in fftool database, the central atom is the third
        '''
        if words[4] not in ['opls']:
            raise Exception('Unsupported improper function: %s' % (words[4]))
        improper = PeriodicImproperTerm(words[2], words[0], words[1], words[3], k=float(words[6]))
        if improper.name in self.improper_terms.keys():
            raise Exception('Duplicated improper term: %s' % self)
        self.improper_terms[improper.name] = improper

    def close(self):
        self._file.close()


class PpfParameterSet(ParameterSet):
    def __init__(self, file, mode='r'):
        super().__init__()
        self.lj_mixing_rule = self.LJ_MIXING_LB
        self.scale_14_vdw = 0.5
        self.scale_14_coulomb = 1.0 / 1.2

        self._file = open(file, mode)
        if mode == 'r':
            self.parse()
        elif mode == 'w':
            raise Exception('Writing is not supported for PpfParameterSet')

    def __del__(self):
        self.close()

    def parse(self):
        '''
        TODO to be implemented
        '''
        pass

    def close(self):
        self._file.close()
