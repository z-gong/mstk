from .ffset import FFSet
from ..ffterm import *
from ..element import Element


class PaduaFFSet(FFSet):
    '''
    fftool of Padua use OPLS convention. Length are in A, energy are in kJ/mol
    The energy term is in k/2 form for bond, angle, dihedral and improper
    '''

    def __init__(self, file):
        super().__init__()
        self.lj_mixing_rule = self.LJ_MIXING_GEOMETRIC
        self.scale_14_vdw = 0.5
        self.scale_14_coulomb = 0.5

        self.parse(file)

    def parse(self, file):
        with open(file) as f:
            lines = f.read().splitlines()
        for line in lines:
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
        atype = AtomType(words[0])
        atype.symbol = Element.guess_from_atom_type(words[0]).symbol
        atype.mass = float(words[2])
        atype.charge = float(words[3])
        (atype.eqt_bond, atype.eqt_ang_c, atype.eqt_ang_s,
         atype.eqt_dih_c, atype.eqt_dih_s,
         atype.eqt_imp_c, atype.eqt_imp_s) = [words[1]] * 7
        if words[4] == 'lj':
            sigma = float(words[5]) / 10  # convert from A to nm
            epsilon = float(words[6])
            lj = LJ126Term(atype.name, atype.name, epsilon, sigma)
            self.vdw_terms[lj.name] = lj
        elif words[4] == 'lj96':
            sigma = float(words[5]) / 10  # convert from A to nm
            epsilon = float(words[6])
            lj = MieTerm(atype.name, atype.name, epsilon, sigma, 9, 6)
            self.vdw_terms[lj.name] = lj
        elif words[4] == 'lj124':
            sigma = float(words[5]) / 10  # convert from A to nm
            epsilon = float(words[6])
            lj = MieTerm(atype.name, atype.name, epsilon, sigma, 12, 4)
            self.vdw_terms[lj.name] = lj
        else:
            raise Exception('Unsupported vdw potential form: %s' % (words[4]))

        if atype.name not in self.atom_types.keys():
            self.atom_types[atype.name] = atype
        else:
            raise Exception('Duplicated atom type: %s' % atype)

    def parse_bond(self, words):
        if words[2] not in ['harm', 'cons']:
            raise Exception('Unsupported bond function: %s' % (words[2]))
        bond = HarmonicBondTerm(words[0], words[1], length=float(words[3]) / 10,
                                k=float(words[4]) * 100 / 2, fixed=(words[2] == 'cons'))
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
        dihedral = OplsDihedralTerm(words[0], words[1], words[2], words[3],
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
                                        phi=180, k=float(words[6]) / 2)
        if improper.name in self.improper_terms.keys():
            raise Exception('Duplicated improper term: %s' % self)
        self.improper_terms[improper.name] = improper

    @staticmethod
    def write(params, file):
        raise Exception('Writing support haven\'t been implemented for fftool format')
