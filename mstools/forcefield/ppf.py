import warnings
from .ffset import FFSet
from .ffterm import *
from .element import Element


class PpfLine():
    def __init__(self, term, key, value, comment):
        self.term = term.strip()
        self.key = key.strip()
        self.value = value.strip()
        self.comment = comment.strip()
        self.name = self.term + '-' + self.key

        version = self.comment.split(',')[0].strip()
        if not version.startswith('V'):
            version = None
        else:
            version = version.lstrip('V').strip()
            try:
                float(version)
            except:
                raise Exception('Invalid version information: %s' % comment)
        self.version = version

    def get_names(self):
        names = [x.strip() for x in self.key.split(',')]
        for i in range(len(names)):
            if names[i] == 'X':
                names[i] = '*'
        return names

    def get_float_values(self):
        return [float(x.strip().strip('*')) for x in self.value.split(',')]


class Ppf(FFSet):
    '''
    In PPF format, there is no 1/2 for all energy terms
    '''

    def __init__(self, *files):
        super().__init__()
        self.lj_mixing_rule = self.LJ_MIXING_LB
        self.vdw_long_range = self.VDW_LONGRANGE_CORRECT
        self.vdw_cutoff = 1.2  # nm
        self.scale_14_vdw = 0.5
        self.scale_14_coulomb = 1.0 / 1.2

        # save all lines to this so that we can keep only the latest version
        self._ppf_lines = {}
        for file in files:
            self.parse(file)

    def parse(self, file):
        with open(file) as f:
            lines = f.read().splitlines()

        for line in lines:
            line = line.strip()
            if line.startswith('#DFF:EQT'):
                section = 'EQT'
                continue
            elif line.startswith('#PROTOCOL'):
                protocol = line.split('=')[-1].strip()
                if protocol != 'AMBER':
                    raise NotImplementedError('Current only AMBER protocol is implemented for PPF')
                continue
            elif line.startswith('#DFF:PPF'):
                section = 'PPF'
                continue
            elif line.startswith('#') or line == '':
                continue

            self._parse_line(line, section)

        self._setup()

    def _parse_line(self, line, section):
        words = [x.strip() for x in line.strip().split(':')]
        if section == 'EQT':
            term = 'EQT'
            key, value = words[:2]
            comment = words[2] if len(words) > 2 else ''
        elif section == 'PPF':
            term, key, value = words[:3]
            comment = words[3] if len(words) > 3 else ''
        else:
            raise Exception('Invalid section: %s' % section)

        ppfLine = PpfLine(term, key, value, comment)
        _ppfLine = self._ppf_lines.get(ppfLine.name)
        if _ppfLine is None:
            self._ppf_lines[ppfLine.name] = ppfLine
        else:
            if ppfLine.version is None or _ppfLine.version is None:
                raise Exception('Duplicated line without version information: %s' % ppfLine.name)
            self._ppf_lines[ppfLine.name] = ppfLine

    def _setup(self):
        for line in self._ppf_lines.values():
            if line.term != 'EQT':
                continue
            atype = AtomType(line.key)
            # _eqt_charge is a temporary variable for matching charge
            (atype.eqt_vdw, atype._eqt_charge, atype.eqt_q_inc, atype.eqt_bond,
             atype.eqt_ang_c, atype.eqt_ang_s, atype.eqt_dih_c, atype.eqt_dih_s,
             atype.eqt_imp_c, atype.eqt_imp_s) = line.value.split()
            atype.version = line.version
            self.atom_types[atype.name] = atype

        for line in self._ppf_lines.values():
            if line.term == 'ATYPE':
                atype = self.atom_types.get(line.key)
                if atype is None:
                    continue
                number, mass = line.get_float_values()
                atype.mass = mass
            elif line.term == 'ATC':
                for atype in self.atom_types.values():
                    if atype._eqt_charge == line.key:
                        atype.charge = line.get_float_values()[0]
            elif line.term == 'BINC':
                term = ChargeIncrementTerm(*line.get_names(), *line.get_float_values())
                term.version = line.version
                self.charge_increment_terms[term.name] = term
            elif line.term == 'N12_6':
                at = line.get_names()[0]
                r_min, epsilon = line.get_float_values()
                term = LJ126Term(at, at, epsilon * 4.184, r_min / 2 ** (1 / 6) / 10)
                term.version = line.version
                self.vdw_terms[term.name] = term
            elif line.term == 'P12_6':
                r_min, epsilon = line.get_float_values()
                term = LJ126Term(*line.get_names(), epsilon * 4.184, r_min / 2 ** (1 / 6) / 10)
                term.version = line.version
                self.pairwise_vdw_terms[term.name] = term
            elif line.term == 'BHARM':
                length, k = line.get_float_values()
                at1, at2 = line.get_names()
                # constrain the bonds involving hydrogen
                fixed = at1.startswith('h_1') or at2.startswith('h_1')
                term = HarmonicBondTerm(at1, at2, length / 10, k * 4.184 * 100, fixed=fixed)
                term.version = line.version
                self.bond_terms[term.name] = term
            elif line.term == 'AHARM':
                theta, k = line.get_float_values()
                at1, at2, at3 = line.get_names()
                # constrain the angles in water
                fixed = at1.startswith('h_1') and at2.startswith('o_2') and at3.startswith('h_1')
                term = HarmonicAngleTerm(at1, at2, at3, theta, k * 4.184, fixed=fixed)
                term.version = line.version
                self.angle_terms[term.name] = term
            elif line.term == 'TCOSP':
                values = line.get_float_values()
                count = len(values) // 3
                term = PeriodicDihedralTerm(*line.get_names())
                for i in range(count):
                    term.add_parameter(
                        values[i * 3], values[i * 3 + 1] * 4.184, int(values[i * 3 + 2]))
                term.version = line.version
                self.dihedral_terms[term.name] = term
            elif line.term == 'IBCOS':
                # the third atom is the center atom for IBCOS improper term
                a1, a2, a3, a4 = line.get_names()
                try:
                    phi, k, n = line.get_float_values()
                except:
                    raise Exception('IBCOS should contain three parameters')
                if phi != 180 or n != 2:
                    raise Exception('IBCOS should have phi=180 and n=2: %s' % line.key)
                term = OplsImproperTerm(a3, a1, a2, a4, k * 4.184)
                term.version = line.version
                self.improper_terms[term.name] = term

    @staticmethod
    def save_to(params: FFSet, file):
        line = '#DFF:EQT\n'
        line += '#AAT :	NB ATC BINC Bond A/C A/S T/C T/S O/C O/S\n'
        for atype in params.atom_types.values():
            line += '%s: %s %s %s %s %s %s %s %s %s %s\n' % (
                atype.name, atype.eqt_vdw, atype.name, atype.eqt_q_inc,
                atype.eqt_bond, atype.eqt_ang_c, atype.eqt_ang_s,
                atype.eqt_dih_c, atype.eqt_dih_s, atype.eqt_imp_c, atype.eqt_imp_s
            )
        line += ('#DFF:PPF\n'
                 '#PROTOCOL = AMBER\n'
                 )
        for atype in params.atom_types.values():
            element = Element.guess_from_atom_type(atype.name)
            line += 'ATYPE: %s: %.5f, %.5f: \n' % (atype.name, element.number, atype.mass)
        for atype in params.atom_types.values():
            line += 'ATC: %s: %.5f: \n' % (atype.name, atype.charge)
        for binc in params.charge_increment_terms.values():
            line += 'BINC: %s, %s: %.5f: ' % (binc.type1, binc.type2, binc.value)
        for vdw in params.vdw_terms.values():
            if isinstance(vdw, LJ126Term):
                line += 'N12_6: %s: %.5f, %.5f: \n' % (
                    vdw.type1, vdw.sigma * 10 * 2 ** (1 / 6), vdw.epsilon / 4.184)
            else:
                raise Exception('Only LJ126Term is implemented')
        for vdw in params.pairwise_vdw_terms.values():
            if isinstance(vdw, LJ126Term):
                line += 'P12_6: %s, %s: %.5f, %.5f: \n' % (
                    vdw.type1, vdw.type2, vdw.sigma * 10 * 2 ** (1 / 6), vdw.epsilon / 4.184)
            else:
                raise Exception('Only LJ126Term is implemented for pairwise vdw terms')
        for bond in params.bond_terms.values():
            if isinstance(bond, HarmonicBondTerm):
                line += 'BHARM: %s, %s: %.5f, %.5f: \n' % (
                    bond.type1, bond.type2, bond.length * 10, bond.k / 4.184 / 100)
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
                line += 'TCOSP: %s, %s, %s, %s: ' % (
                    dihedral.type1, dihedral.type2, dihedral.type3, dihedral.type4)
                str_paras = []
                for par in dihedral.parameters:
                    str_paras.append('%.1f, %.5f, %i' % (par.phi, par.k / 4.184, par.n))
                line += ', '.join(str_paras) + ': \n'
            else:
                raise Exception('Only PeriodicDihedralTerm is implemented')
        for improper in params.improper_terms.values():
            if isinstance(improper, OplsImproperTerm):
                # the third atom is the center atom for IBCOS improper term
                line += 'IBCOS: %s, %s, %s, %s: 180, %.5f, 2: \n' % (
                    improper.type2, improper.type3, improper.type1, improper.type4,
                    improper.k / 4.184)
            else:
                raise Exception('Only PeriodicImproperTerm is implemented')

        if len(params.polarizable_terms) > 0:
            warnings.warn('Polarizable parameters are ignored because PPF does not support them')

        with open(file, 'w') as f:
            f.write(line)
