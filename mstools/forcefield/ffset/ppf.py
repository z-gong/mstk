from .ffset import FFSet
from ..ffterm import *
from ..element import Element


class PpfLine():
    def __init__(self, term, key, value, comment):
        self.term = term.strip()
        self.key = key.strip()
        self.value = value.strip()
        self.comment = comment.strip()
        self.name = self.term + '-' + self.key

        try:
            version = comment.split(',')[0].strip('V').strip()
        except:
            version = None
        if version is not None:
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


class PpfFFSet(FFSet):
    '''
    In PPF format, there is no 1/2 for all energy terms
    '''

    def __init__(self, file):
        super().__init__()
        self.lj_mixing_rule = self.LJ_MIXING_LB
        self.scale_14_vdw = 0.5
        self.scale_14_coulomb = 1.0 / 1.2

        self.parse(file)

    def parse(self, ppf_file):
        # save all lines to this so that we can keep only the latest version
        self._ppf_lines = {}

        with open(ppf_file) as f:
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
            (atype.eqt_vdw, atype._eqt_charge, atype.eqt_charge_increment, atype.eqt_bond,
             atype.eqt_angle_center, atype.eqt_angle_side, atype.eqt_dihedral_center,
             atype.eqt_dihedral_side, atype.eqt_improper_center,
             atype.eqt_improper_side) = line.value.split()
            atype.version = line.version
            self.atom_types[atype.name] = atype

        for line in self._ppf_lines.values():
            if line.term == 'ATYPE':
                atype = self.atom_types.get(line.key)
                if atype is None:
                    continue
                number, mass = line.get_float_values()
                atype.symbol = Element(int(number)).symbol
                atype.mass = mass
            elif line.term == 'ATC':
                for atype in self.atom_types.values():
                    if atype._eqt_charge == line.key:
                        atype.charge = line.get_float_values()[0]
                        atype._eqt_charge = atype.name
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
                term = HarmonicBondTerm(*line.get_names(), length / 10, k * 4.184 * 100)
                term.version = line.version
                self.bond_terms[term.name] = term
            elif line.term == 'AHARM':
                theta, k = line.get_float_values()
                term = HarmonicAngleTerm(*line.get_names(), theta, k * 4.184)
                term.version = line.version
                self.angle_terms[term.name] = term
            elif line.term == 'TCOSP':
                values = line.get_float_values()
                phi_list = values[0::3]
                k_list = values[1::3]
                n_list = list(map(int, values[2::3]))
                term = PeriodicDihedralTerm(*line.get_names(), n_list, k_list, phi_list)
                term.version = line.version
                self.dihedral_terms[term.name] = term
            elif line.term == 'IBCOS':
                # the third atom is the center atom in PPF
                a1, a2, a3, a4 = line.get_names()
                phi, k, multiplicity = line.get_float_values()
                if multiplicity != 2:
                    raise Exception('Multiplicity in IBCOS should equal to 2: %s' % line.key)
                term = PeriodicImproperTerm(a3, a1, a2, a4, phi, k * 4.184)
                term.version = line.version
                self.improper_terms[term.name] = term

    @staticmethod
    def write(params: FFSet, file):
        line = '#DFF:EQT\n'
        line += '#AAT :	NB ATC BINC Bond A/C A/S T/C T/S O/C O/S\n'
        for atype in params.atom_types.values():
            line += '%s: %s %s %s %s %s %s %s %s %s %s\n' % (
                atype.name, atype.eqt_vdw, atype.name, atype.eqt_charge,
                atype.eqt_bond, atype.eqt_angle_center, atype.eqt_angle_side,
                atype.eqt_dihedral_center, atype.eqt_dihedral_side,
                atype.eqt_improper_center, atype.eqt_improper_side
            )
        line += ('#DFF:PPF\n'
                 '#PROTOCOL = AMBER\n'
                 )
        for atype in params.atom_types.values():
            element = Element(atype.symbol)
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
                    str_paras.append('%.1f, %.5f, %i' % (par.phi, par.k / 4.184, par.multiplicity))
                line += ', '.join(str_paras) + ': \n'
            elif isinstance(dihedral, OplsDihedralTerm):
                line += 'TCOSP: %s, %s, %s, %s: ' % (
                    dihedral.type1, dihedral.type2, dihedral.type3, dihedral.type4)
                str_paras = []
                if dihedral.k1 != 0: str_paras.append('%.1f, %.5f, 1' % (0.0, dihedral.k1 / 4.184))
                if dihedral.k2 != 0: str_paras.append('%.1f, %.5f, 2' % (180, dihedral.k2 / 4.184))
                if dihedral.k3 != 0: str_paras.append('%.1f, %.5f, 3' % (0.0, dihedral.k3 / 4.184))
                if dihedral.k4 != 0: str_paras.append('%.1f, %.5f, 4' % (180, dihedral.k4 / 4.184))
                if {dihedral.k1, dihedral.k2, dihedral.k3, dihedral.k4} == {0}:
                    str_paras.append('0, 0 , 1')
                line += ', '.join(str_paras) + ': \n'
            else:
                raise Exception('Only PeriodicDihedralTerm and OplsDihedralTerm are implemented')
        for improper in params.improper_terms.values():
            if isinstance(improper, PeriodicImproperTerm):
                line += 'IBCOS: %s, %s, %s, %s: %.1f, %.5f, 2: \n' % (
                    improper.type2, improper.type3, improper.type1, improper.type4,
                    improper.phi, improper.k / 4.184)
            else:
                raise Exception('Only PeriodicImproperTerm is implemented')

        with open(file, 'w') as f:
            f.write(line)
