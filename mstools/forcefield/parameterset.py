import math
from .ffterm import *
from .element import Element


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
        atom_types = list(self.atom_types.values())
        for i in range(len(atom_types)):
            for j in range(i, len(atom_types)):
                vdw = self.get_vdw_term(atom_types[i], atom_types[j])
                self.vdw_terms[vdw.name] = vdw

    def get_vdw_term(self, type1: AtomType, type2: AtomType):
        '''
        get vdW term between two atom types
        if not exist, assuming it's LJ126 and generate it using combination rule
        '''
        atype1 = type1.eqt_vdw
        atype2 = type2.eqt_vdw
        vdw = self.vdw_terms.get(VdwTerm(atype1.name, atype2.name).name)
        if vdw is not None:
            return vdw

        lj = LJ126Term(atype1.name, atype2.name, 0.0, 1.0)
        if self.lj_mixing_rule == self.LJ_MIXING_LB:
            lj.epsilon = math.sqrt(atype1.epsilon * atype2.epsilon)
            lj.sigma = (atype1.sigma + atype2.sigma) / 2
        elif self.lj_mixing_rule == self.LJ_MIXING_GEOMETRIC:
            lj.epsilon = math.sqrt(atype1.epsilon * atype2.epsilon)
            lj.sigma = math.sqrt(atype1.sigma * atype2.sigma)
        else:
            raise Exception('Unknown mixing rule for LJ126 parameters')

        return lj


class FFToolParameterSet(ParameterSet):
    '''
    FFTool use OPLS convention. Length are in A, energy are in kJ/mol
    The energy term is in k/2 form for bond, angle, dihedral and improper
    '''

    def __init__(self, file):
        super().__init__()
        self.lj_mixing_rule = self.LJ_MIXING_GEOMETRIC
        self.scale_14_vdw = 0.5
        self.scale_14_coulomb = 0.5

        self._file = open(file)
        self.parse()
        self._file.close()

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

        self.setup_eqts()

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

    def setup_eqts(self):
        for name, atype in self._eq_atom_types.items():
            if name not in self.atom_types:
                self.atom_types[name] = atype

        for atype in self.atom_types.values():
            atype.eqt_bond = self.atom_types[atype._eq_type_name]
            atype.eqt_angle_center = self.atom_types[atype._eq_type_name]
            atype.eqt_angle_side = self.atom_types[atype._eq_type_name]
            atype.eqt_dihedral_center = self.atom_types[atype._eq_type_name]
            atype.eqt_dihedral_side = self.atom_types[atype._eq_type_name]
            atype.eqt_improper_center = self.atom_types[atype._eq_type_name]
            atype.eqt_improper_side = self.atom_types[atype._eq_type_name]

    def parse_bond(self, words):
        if words[2] not in ['harm', 'cons']:
            raise Exception('Unsupported bond function: %s' % (words[2]))
        bond = HarmonicBondTerm(words[0], words[1], length=float(words[3]) / 10,
                                k=float(words[4]) / 2, fixed=(words[2] == 'cons'))
        if bond.name in self.bond_terms.keys():
            raise Exception('Duplicated bond term: %s' % self)
        self.bond_terms[bond.name] = bond

    def parse_angle(self, words):
        if words[3] not in ['harm']:
            raise Exception('Unsupported angle function: %s' % (words[3]))
        angle = HarmonicAngleTerm(words[0], words[1], words[2], theta=float(words[4]),
                                  k=float(words[5]) / 2)
        if angle.name in self.angle_terms.keys():
            raise Exception('Duplicated angle term: %s' % self)
        self.angle_terms[angle.name] = angle

    def parse_dihedral(self, words):
        if words[4] not in ['opls']:
            raise Exception('Unsupported dihedral function: %s' % (words[4]))
        dihedral = FourierDihedralTerm(words[0], words[1], words[2], words[3],
                                       k1=float(words[5]) / 2, k2=float(words[6]) / 2,
                                       k3=float(words[7]) / 2, k4=float(words[8]) / 2
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
        improper = PeriodicImproperTerm(words[2], words[0], words[1], words[3],
                                        k=float(words[6]) / 2)
        if improper.name in self.improper_terms.keys():
            raise Exception('Duplicated improper term: %s' % self)
        self.improper_terms[improper.name] = improper

    def write(self, file):
        '''
        Writing is not supported for FFTool
        '''


class PpfParameterSet(ParameterSet):
    '''
    In pff format, there is no 1/2 for all energy terms
    '''

    def __init__(self, file):
        super().__init__()
        self.lj_mixing_rule = self.LJ_MIXING_LB
        self.scale_14_vdw = 0.5
        self.scale_14_coulomb = 1.0 / 1.2

        self._file = open(file)
        self.parse()
        self._file.close()

    def parse(self):
        '''
        TODO to be implemented
        '''
        pass

    @staticmethod
    def write(params, file):
        line = '#DFF:EQT\n'
        line += '#AAT :	NB ATC BINC Bond A/C A/S T/C T/S O/C O/S\n'
        for atype in params.atom_types.values():
            line += '%s: %s %s %s %s %s %s %s %s %s %s\n' % (
                atype.name, atype.eqt_vdw.name, atype.eqt_charge.name, atype.eqt_charge.name,
                atype.eqt_bond.name,
                atype.eqt_angle_center.name, atype.eqt_angle_side.name,
                atype.eqt_dihedral_center.name, atype.eqt_dihedral_side.name,
                atype.eqt_improper_center.name, atype.eqt_improper_side.name
            )
        line += ('#DFF:PPF\n'
                 '#PROTOCOL = AMBER\n'
                 '#TYPINGRULE = \n'
                 )
        for atype in params.atom_types.values():
            element = Element.guess_from_atom_type(atype.name)
            line += 'ATYPE: %s: %.5f, %.5f: \n' % (atype.name, element.number, element.mass)
        for atype in params.atom_types.values():
            line += 'ATC: %s: %.5f: \n' % (atype.name, atype.charge)
        for atype in params.atom_types.values():
            line += 'N12_6: %s: %.5f, %.5f: \n' % (
                atype.name, atype.sigma * 10 * 2 ** (1 / 6), atype.epsilon / 4.184)
        for binc in params.charge_increment_terms.values():
            line += 'BINC: %s, %s: %.5f: ' % (binc.type1, binc.type2, binc.value)
        for bond in params.bond_terms.values():
            if isinstance(bond, HarmonicBondTerm):
                line += 'BHARM: %s, %s: %.5f, %.5f: \n' % (
                    bond.type1, bond.type2, bond.length * 10, bond.k / 4.184)
            else:
                raise Exception('Only HarmonicBondTerm is implemented')
        for angle in params.angle_terms.values():
            if isinstance(angle, HarmonicAngleTerm):
                line += 'AHARM: %s, %s, %s: %.5f, %.5f: \n' % (
                    angle.type1, angle.type2, angle.type3, angle.theta, angle.k / 4.184)
            else:
                raise Exception('Only HarmonicAngleTerm is implemented')
        for dihedral in params.dihedral_terms.values():
            if isinstance(dihedral, PeriodicDihedralTerm):
                line += 'TCOSP: %s, %s, %s, %s:' % (
                    dihedral.type1, dihedral.type2, dihedral.type3, dihedral.type4)
                for par in dihedral.parameters:
                    line += f' {par.phi}, {par.k}, {par.multiplicity}'
                line += ' : \n'
            elif isinstance(dihedral, FourierDihedralTerm):
                # if dihedral.k4 != 0:
                #     print(f'4th terms of {str(dihedral)} is ignored')
                line += 'TCOSP: %s, %s, %s, %s: ' % (
                    dihedral.type1, dihedral.type2, dihedral.type3, dihedral.type4)
                str_paras = []
                if dihedral.k1 != 0: str_paras.append('0.0, %.5f, 1' % (dihedral.k1 / 4.184))
                if dihedral.k2 != 0: str_paras.append('180, %.5f, 2' % (dihedral.k2 / 4.184))
                if dihedral.k3 != 0: str_paras.append('0.0, %.5f, 3' % (dihedral.k3 / 4.184))
                if dihedral.k4 != 0: str_paras.append('180, %.5f, 4' % (dihedral.k4 / 4.184))
                if {dihedral.k1, dihedral.k2, dihedral.k3, dihedral.k4} == {0}:
                    str_paras.append('0.0, 0 , 1')
                line += ', '.join(str_paras) + ': \n'
            else:
                raise Exception('Only PeriodicDihedralTerm is implemented')
        for improper in params.improper_terms.values():
            if isinstance(improper, PeriodicImproperTerm):
                line += 'IBCOS: %s, %s, %s, %s: 180, %.5f, 2: \n' % (
                    improper.type2, improper.type3, improper.type1, improper.type4,
                    improper.k / 4.184)
            else:
                raise Exception('Only PeriodicImproperTerm is implemented')

        with open(file, 'w') as f:
            f.write(line)
