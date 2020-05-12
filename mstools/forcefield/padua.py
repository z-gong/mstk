import itertools
from .ffset import FFSet
from .ffterm import *
from .element import Element
from .. import logger


class Padua(FFSet):
    '''
    fftool of Padua use OPLS convention. Length are in A, energy are in kJ/mol
    The energy term is in k/2 form for bond, angle, dihedral and improper
    '''

    def __init__(self, *files):
        super().__init__()
        self.lj_mixing_rule = self.LJ_MIXING_GEOMETRIC
        self.vdw_long_range = self.VDW_LONGRANGE_CORRECT
        self.vdw_cutoff = 1.2  # nm
        self.scale_14_vdw = 0.5
        self.scale_14_coulomb = 0.5

        self.ljscaler: PaduaLJScaler = None

        for file in files:
            self._parse(file)

        if self.ljscaler is not None:
            self.ljscaler.scale(self)

    def _parse(self, file):
        with open(file) as f:
            lines = f.read().splitlines()

        _section = ''
        for line in lines:
            if line.startswith('#') or line.strip() == '':
                continue

            if line.lower().startswith('atom'):
                _section = 'atoms'
                continue
            elif line.lower().startswith('bond'):
                _section = 'bonds'
                continue
            elif line.lower().startswith('angl'):
                _section = 'angles'
                continue
            elif line.lower().startswith('dihe'):
                _section = 'dihedrals'
                continue
            elif line.lower().startswith('impro'):
                _section = 'impropers'
                continue
            elif line.lower().startswith('polar'):
                _section = 'polarizations'
                continue
            elif line.lower().startswith('scale_sigma') or line.lower().startswith('monomer'):
                self.ljscaler = PaduaLJScaler(file)
                return

            words = line.strip().split()
            if _section == 'atoms':
                self._parse_atom(words)
            elif _section == 'bonds':
                self._parse_bond(words)
            elif _section == 'angles':
                self._parse_angle(words)
            elif _section == 'dihedrals':
                self._parse_dihedral(words)
            elif _section == 'impropers':
                self._parse_improper(words)
            elif _section == 'polarizations':
                self._parse_polarization(words)
            else:
                raise Exception('Invalid header for fftool: %s' % line)

        # TODO not robust enough
        # If this is a polarizable FF, add an AtomType and VdwTerm for Drude particles
        if len(self.polarizable_terms) > 0:
            type_drude = 'DP_'
            if type_drude not in self.atom_types:
                dtype = AtomType(type_drude)
                vdw = LJ126Term(type_drude, type_drude, 0.0, 0.0)
                self.add_term(dtype)
                self.add_term(vdw)

        # If H* is defined in polarizable term, alpha of H will be merged into attached heavy atoms and the term will be removed
        if 'H*' in self.polarizable_terms:
            hterm = self.polarizable_terms.pop('H*')
            for pterm in self.polarizable_terms.values():
                pterm.merge_alpha_H = hterm.alpha
            logger.warning(f'{str(hterm)} found in polarizable term. '
                           f'Polarizability of H will be merged into attached heavy atoms')

    def _parse_atom(self, words):
        atype = AtomType(words[0])
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
        else:
            raise Exception('Unsupported vdw potential form: %s' % (words[4]))

        if atype.name not in self.atom_types.keys():
            self.atom_types[atype.name] = atype
        else:
            raise Exception('Duplicated atom type: %s' % str(atype))

    def _parse_bond(self, words):
        if words[2] not in ['harm', 'cons']:
            raise Exception('Unsupported bond function: %s' % (words[2]))
        term = HarmonicBondTerm(words[0], words[1], length=float(words[3]) / 10,
                                k=float(words[4]) * 100 / 2, fixed=(words[2] == 'cons'))
        if term.name in self.bond_terms.keys():
            raise Exception('Duplicated bond term: %s' % str(term))
        self.bond_terms[term.name] = term

    def _parse_angle(self, words):
        if words[3] not in ['harm', 'cons']:
            raise Exception('Unsupported angle function: %s' % (words[3]))
        term = HarmonicAngleTerm(words[0], words[1], words[2], theta=float(words[4]),
                                 k=float(words[5]) / 2, fixed=(words[3] == 'cons'))
        if term.name in self.angle_terms.keys():
            raise Exception('Duplicated angle term: %s' % str(term))
        self.angle_terms[term.name] = term

    def _parse_dihedral(self, words):
        if words[4] not in ['opls']:
            raise Exception('Unsupported dihedral function: %s' % (words[4]))
        term = PeriodicDihedralTerm(words[0], words[1], words[2], words[3])
        k1, k2, k3, k4 = list(map(lambda x: float(x) / 2, words[5:9]))
        if k1 != 0: term.add_parameter(0.0, k1, 1)
        if k2 != 0: term.add_parameter(180, k2, 2)
        if k3 != 0: term.add_parameter(0.0, k3, 3)
        if k4 != 0: term.add_parameter(180, k4, 4)
        if term.name in self.dihedral_terms.keys():
            raise Exception('Duplicated dihedral term: %s' % str(term))
        self.dihedral_terms[term.name] = term

    def _parse_improper(self, words):
        '''
        For improper i-j-k-l
        mstools always use CHARMM convention, i is the central atom
        fftool use OPLS convention, k is the central atom
        '''
        if words[4] not in ['opls']:
            raise Exception('Unsupported improper function: %s' % (words[4]))
        if {float(words[5]), float(words[7]), float(words[8])} != {0}:
            raise Exception('k1, k3 and k4 for opls improper should equal to 0')
        term = OplsImproperTerm(words[2], words[0], words[1], words[3], k=float(words[6]) / 2)
        if term.name in self.improper_terms.keys():
            raise Exception('Duplicated improper term: %s' % str(term))
        self.improper_terms[term.name] = term

    def _parse_polarization(self, words):
        name, mass, charge, k, alpha, thole = words
        term = DrudeTerm(name, float(alpha) / 1000, float(thole))
        term.mass = float(mass)
        term.k = float(k) / 2 * 100  # convert from kJ/mol/A^2 to kJ/mol/nm^2
        if term.name in self.polarizable_terms.keys():
            raise Exception('Duplicated drude term: %s' % str(term))
        self.polarizable_terms[term.name] = term


class Monomer():
    def __init__(self, name: str, charge: float, dipole: float, alpha: float = None):
        self.name = name
        self.charge = charge
        self.dipole = dipole
        self.alpha = alpha or -1  # set to -1 if not provided
        self.atoms: [str] = []  # list of atom types

    def __repr__(self):
        return '<Monomer: %s>' % self.name


class Dimer():
    C0 = 0.25
    C1 = 0.11

    def __init__(self, m1: Monomer, m2: Monomer, distance: float, scale_factor: float = None):
        '''
        If scale factor is None, then need to calculate scale factor later
        based on charge, dipole, polarizability and distance of monomers
        '''
        self.monomer1 = m1
        self.monomer2 = m2
        self.distance = distance
        self.scale_factor = scale_factor

    def __repr__(self):
        return '<Dimer: %s %s>' % (self.monomer1.name, self.monomer2.name)

    def predict_scale_epsilon(self):
        m1, m2 = self.monomer1, self.monomer2
        if (m1.charge != 0 or m1.dipole != 0) and m1.alpha == -1:
            raise Exception(
                'Error in predicting scaling factor: alpha required for charged or polar fragment: %s' % m1.name)
        if (m2.charge != 0 or m2.dipole != 0) and m2.alpha == -1:
            raise Exception(
                'Error in predicting scaling factor: alpha required for charged or polar fragment: %s' % m2.name)

        k = self.C0 * self.distance ** 2 * (m1.charge ** 2 / m1.alpha + m2.charge ** 2 / m2.alpha) + \
            self.C1 * (m1.dipole ** 2 / m1.alpha + m2.dipole ** 2 / m2.alpha)
        self.scale_factor = round(1 / (1 + k), 3)


class PaduaLJScaler():
    def __init__(self, file):
        self.monomers = []
        self.dimers = []
        self.scale_sigma = 1.0
        self._file = file
        self._parse(file)

    def __repr__(self):
        return '<PaduaLJScaler: %s>' % self._file

    def scale(self, ffset: FFSet):
        # must scale pairwise vdW terms first
        # because they are generated from self vdW terms by combination rule
        # the scaled pairwise terms must be added into the ff set
        for type1, type2 in itertools.combinations(ffset.atom_types.values(), 2):
            vdw = ffset.get_vdw_term(type1, type2)
            if type(vdw) == LJ126Term:
                self.scale_vdw(vdw)
                if vdw.name not in ffset.pairwise_vdw_terms:
                    ffset.add_term(vdw)
        for vdw in ffset.vdw_terms.values():
            if type(vdw) == LJ126Term:
                self.scale_vdw(vdw)

    def _parse(self, file):
        '''
        Scaling factor can be provided in DIMERS seciton
        K_SAPT2 section is ignored
        :param file:
        :return:
        '''
        with open(file) as f:
            lines = f.read().splitlines()

        _section = ''
        for line in lines:
            if line.strip() == '' or line.startswith('#'):
                continue

            words = line.strip().split()
            if words[0] == 'MONOMERS':
                _section = 'MONOMERS'
                continue
            elif words[0] == 'DIMERS':
                _section = 'DIMERS'
                continue
            elif words[0] == 'K_SAPT2':
                _section = 'K_SAPT2'
                continue
            elif words[0] == 'ATOMS':
                _section = 'ATOMS'
                continue
            elif words[0] == 'SCALE_SIGMA':
                self.scale_sigma = float(words[1])
                continue

            if _section == 'MONOMERS':
                name = words[0].strip('+').strip('-')
                alpha = None
                if len(words) > 3:
                    alpha = float(words[3])
                self.monomers.append(Monomer(name, float(words[1]), float(words[2]), alpha))

            elif _section == 'DIMERS':
                m1 = next(m for m in self.monomers if m.name == words[0])
                m2 = next(m for m in self.monomers if m.name == words[1])
                scale = None
                if len(words) > 3:
                    scale = float(words[3])
                self.dimers.append(Dimer(m1, m2, float(words[2]), scale))

            elif _section == 'ATOMS':
                monomer = next(m for m in self.monomers if m.name == words[0])
                monomer.atoms = words[1:]

    def predict_scale_epsilon(self, atom1, atom2):
        '''
        :param atom1:
        :param atom2:
        :return: k = None means monomer or dimer data not found. will not scale
                 k = positive float means scale normally
        '''
        m1 = m2 = None
        for monomer in self.monomers:
            if atom1 in monomer.atoms:
                m1 = monomer
            if atom2 in monomer.atoms:
                m2 = monomer
        if m1 is None or m2 is None:
            return None

        for dimer in self.dimers:
            if (m1 == dimer.monomer1 and m2 == dimer.monomer2) or (
                    m2 == dimer.monomer1 and m1 == dimer.monomer2):
                dimer = dimer
                break
        else:
            return None

        if dimer.scale_factor is None:
            dimer.predict_scale_epsilon()

        return dimer.scale_factor

    def scale_vdw(self, term: LJ126Term):
        k = self.predict_scale_epsilon(term.type1, term.type2)
        if k is None:
            return

        term.epsilon *= k
        term.comments.append('eps*%.3f' % k)

        factor = self.scale_sigma
        if factor == 1.0:
            return
        term.sigma *= factor
        term.comments.append('sig*%.3f' % factor)
