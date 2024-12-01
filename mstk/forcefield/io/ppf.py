from mstk.forcefield.forcefield import ForceField
from mstk.forcefield.ffterm import *
from mstk.chem.element import Element
from mstk.forcefield.errors import *
from mstk import logger


class PpfLine:
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


class Ppf:
    '''
    Parse ForceField from PPF file.

    PPF is the native format of `Direct Force Field`.
    Currently, only the DFF-AMBER functional form of PPF is supported,
    which contains following terms:

    * :class:`~mstk.forcefield.AtomType` with equivalent tables, masses and partial charges
    * :class:`~mstk.forcefield.ChargeIncrementTerm`
    * :class:`~mstk.forcefield.HarmonicBondTerm`
    * :class:`~mstk.forcefield.HarmonicAngleTerm`
    * :class:`~mstk.forcefield.PeriodicDihedralTerm`
    * :class:`~mstk.forcefield.PeriodicImproperTerm`

    In PPF file, there is no leading 1/2 for harmonic energy terms.
    The length is in unit of Angstrom, and angle is in unit of degree.
    The energies for bond, angle and vdW are in unit of kcal/mol/A^2, kcal/mol/rad^2 and kcal/mol.
    The two values for LJ interaction are r_min and epsilon.

    PPF format can optionally store the version of parameters.
    If a force field term describing one topology element appear several times,
    only the one with the latest version will be loaded.
    If there are duplicated terms without version record, an Exception will be raised.

    Parameters
    ----------
    file : str

    Attributes
    ----------
    forcefield : ForceField
    '''

    def __init__(self, file):
        ff = ForceField()
        ff.lj_mixing_rule = ForceField.LJ_MIXING_LB
        ff.vdw_long_range = ForceField.VDW_LONGRANGE_CORRECT
        ff.vdw_cutoff = 1.2  # nm
        ff.scale_14_vdw = 0.5
        ff.scale_14_coulomb = 1.0 / 1.2

        # save all lines to this so that we can keep only the latest version
        self._ppf_lines = {}
        self._parse(file)
        self._setup(ff)

        self.forcefield = ff

    def _parse(self, file):
        with open(file) as f:
            lines = f.read().splitlines()

        _section = ''
        for line in lines:
            line = line.strip()
            if line.startswith('#DFF:EQT'):
                _section = 'EQT'
                continue
            elif line.startswith('#PROTOCOL'):
                protocol = line.split('=')[-1].strip()
                if protocol != 'AMBER':
                    raise NotImplementedError('Current only AMBER protocol is implemented for PPF')
                continue
            elif line.startswith('#DFF:PPF'):
                _section = 'PPF'
                continue
            elif line.startswith('#') or line == '':
                continue
            if _section == '':
                raise Exception('Invalid PPF file')

            self._parse_line(line, _section)

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

        ppf_line = PpfLine(term, key, value, comment)
        # only keep the line of the latest version
        old_line = self._ppf_lines.get(ppf_line.name)
        if old_line is None:
            self._ppf_lines[ppf_line.name] = ppf_line
        else:
            if ppf_line.version is None or old_line.version is None:
                raise Exception('Duplicated line without version information: %s' % ppf_line.name)
            if float(ppf_line.version) > float(old_line.version):
                self._ppf_lines[ppf_line.name] = ppf_line

    def _setup(self, ff):
        for line in self._ppf_lines.values():
            if line.term != 'EQT':
                continue
            atype = AtomType(line.key)
            # _eqt_charge_ is a temporary variable for matching charge
            (atype.eqt_vdw, atype._eqt_charge_, atype.eqt_bci, atype.eqt_bond,
             atype.eqt_ang_c, atype.eqt_ang_s, atype.eqt_dih_c, atype.eqt_dih_s,
             atype.eqt_imp_c, atype.eqt_imp_s) = line.value.split()
            atype.version = line.version
            ff.atom_types[atype.name] = atype

        for line in self._ppf_lines.values():
            if line.term == 'ATYPE':
                atype = ff.atom_types.get(line.key)
                if atype is None:
                    continue
                number, mass = line.get_float_values()
                atype.mass = mass
            elif line.term == 'ATC':
                for atype in ff.atom_types.values():
                    if atype._eqt_charge_ == line.key:
                        atype.charge = line.get_float_values()[0]
            elif line.term == 'BINC':
                try:
                    term = ChargeIncrementTerm(*line.get_names(), *line.get_float_values())
                except ChargeIncrementNonZeroError:
                    logger.warning('Invalid BINC line ignored: %s' % line.key)
                else:
                    term.version = line.version
                    ff.bci_terms[term.name] = term
            elif line.term == 'N12_6':
                at = line.get_names()[0]
                r_min, epsilon = line.get_float_values()
                term = LJ126Term(at, at, epsilon * 4.184, r_min / 2 ** (1 / 6) / 10)
                term.version = line.version
                ff.vdw_terms[term.name] = term
            elif line.term == 'P12_6':
                r_min, epsilon = line.get_float_values()
                term = LJ126Term(*line.get_names(), epsilon * 4.184, r_min / 2 ** (1 / 6) / 10)
                term.version = line.version
                ff.pairwise_vdw_terms[term.name] = term
            elif line.term == 'BHARM':
                length, k = line.get_float_values()
                at1, at2 = line.get_names()
                # constrain the bonds involving hydrogen
                fixed = at1.startswith('h_1') or at2.startswith('h_1')
                term = HarmonicBondTerm(at1, at2, length / 10, k * 4.184 * 100, fixed=fixed)
                term.version = line.version
                ff.bond_terms[term.name] = term
            elif line.term == 'AHARM':
                theta, k = line.get_float_values()
                at1, at2, at3 = line.get_names()
                # constrain the angles in water
                fixed = at1.startswith('h_1') and at2.startswith('o_2') and at3.startswith('h_1')
                term = HarmonicAngleTerm(at1, at2, at3, theta * DEG2RAD, k * 4.184, fixed=fixed)
                term.version = line.version
                ff.angle_terms[term.name] = term
            elif line.term == 'TCOSP':
                values = line.get_float_values()
                count = len(values) // 3
                term = PeriodicDihedralTerm(*line.get_names())
                for i in range(count):
                    phi, k, n = values[i * 3], values[i * 3 + 1] * 4.184, int(values[i * 3 + 2])
                    # ignore nonsense values
                    if n == 0 and k == 0:
                        continue
                    term.add_phase(phi * DEG2RAD, k, n)
                term.version = line.version
                ff.dihedral_terms[term.name] = term
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
                ff.improper_terms[term.name] = term

        # update the mass based on eqt_vdw
        for atype in ff.atom_types.values():
            if atype.eqt_vdw != atype.name:
                atype_eq = ff.atom_types.get(atype.eqt_vdw)
                if atype_eq is None:
                    raise Exception('vdW equivalent type for %s not found in PPF' % (str(atype)))
                atype.mass = atype_eq.mass


ForceField.register_format('.ppf', Ppf)
