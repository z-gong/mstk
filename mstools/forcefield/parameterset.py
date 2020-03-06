from .ffterm import *


class ParameterSet():
    def __init__(self):
        self.atom_types = {}
        self.bond_terms = {}
        self.angle_terms = {}
        self.dihedral_terms = {}
        self.improper_terms = {}
        self.nonbonded_terms = {}


class FFToolParameterSet(ParameterSet):
    def __init__(self, file, mode='r'):
        super().__init__()
        self._file = open(file, mode)
        if mode == 'r':
            self.parse()
        elif mode == 'w':
            raise Exception('Writing is not supported for FFToolParameterSet')

    def __del__(self):
        self.close()

    def parse(self):
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

    def parse_atom(self, words):
        atype_eq = AtomType(words[1])
        if atype_eq.key not in self.atom_types.keys():
            self.atom_types[atype_eq.key] = atype_eq

        atype = AtomType(words[0], mass=float(words[2]))
        atype.charge = float(words[3])
        if words[4] == 'lj':
            atype.sigma = float(words[5]) / 10  # convert from A to nm
            atype.epsilon = float(words[6])
        else:
            raise Exception('Unsupported nonbonded function: %s' % (words[4]))

        if atype.key in self.atom_types.keys() and atype.key != atype_eq.key:
            raise Exception('Duplicated atom type: %s' % atype)
        self.atom_types[atype.key] = atype

        atype.equivalent_type_bond = self.atom_types[atype_eq.key]
        atype.equivalent_type_angle_center = self.atom_types[atype_eq.key]
        atype.equivalent_type_angle_side = self.atom_types[atype_eq.key]
        atype.equivalent_type_dihedral_center = self.atom_types[atype_eq.key]
        atype.equivalent_type_dihedral_side = self.atom_types[atype_eq.key]
        atype.equivalent_type_improper_center = self.atom_types[atype_eq.key]
        atype.equivalent_type_improper_side = self.atom_types[atype_eq.key]

    def parse_bond(self, words):
        if words[2] not in ['harm', 'cons']:
            raise Exception('Unsupported bond function: %s' % (words[2]))
        bond = HarmonicBondTerm(words[0], words[1], length=float(words[3]) / 10, k=float(words[4]))
        if bond.key in self.bond_terms.keys():
            raise Exception('Duplicated bond term: %s' % self)
        self.bond_terms[bond.key] = bond

    def parse_angle(self, words):
        if words[3] not in ['harm']:
            raise Exception('Unsupported angle function: %s' % (words[3]))
        angle = HarmonicAngleTerm(words[0], words[1], words[2], theta=float(words[4]), k=float(words[5]))
        if angle.key in self.angle_terms.keys():
            raise Exception('Duplicated angle term: %s' % self)
        self.angle_terms[angle.key] = angle

    def parse_dihedral(self, words):
        if words[4] not in ['opls']:
            raise Exception('Unsupported dihedral function: %s' % (words[4]))
        dihedral = FourierDihedralTerm(words[0], words[1], words[2], words[3],
                                       k1=float(words[5]), k2=float(words[6]),
                                       k3=float(words[7]), k4=float(words[8])
                                       )
        if dihedral.key in self.dihedral_terms.keys():
            raise Exception('Duplicated dihedral term: %s' % self)
        self.dihedral_terms[dihedral.key] = dihedral

    def parse_improper(self, words):
        '''
        in fftool database, the central atom is the third
        '''
        if words[4] not in ['opls']:
            raise Exception('Unsupported improper function: %s' % (words[4]))
        improper = PeriodicImproperTerm(words[2], words[0], words[1], words[3], k=float(words[6]))
        if improper.key in self.improper_terms.keys():
            raise Exception('Duplicated improper term: %s' % self)
        self.improper_terms[improper.key] = improper

    def close(self):
        self._file.close()


class PpfParameterSet(ParameterSet):
    def __init__(self, file, mode='r'):
        super().__init__()
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
