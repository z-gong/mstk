import itertools
from mstk.forcefield.forcefield import ForceField
from mstk.forcefield.ffterm import *
from mstk.forcefield.errors import *
from mstk import logger


class Padua:
    '''
    Generate ForceField from the force field file of `fftool`.

    The `fftool` developed by Agilio Pauda use OPLS convention with extension of Drude polarization,
    which contains following terms:

    * :class:`~mstk.forcefield.AtomType` with equivalent tables, masses and partial charges
    * :class:`~mstk.forcefield.HarmonicBondTerm`
    * :class:`~mstk.forcefield.HarmonicAngleTerm`
    * :class:`~mstk.forcefield.PeriodicDihedralTerm`
    * :class:`~mstk.forcefield.OplsImproperTerm`
    * :class:`~mstk.forcefield.DrudePolarTerm`

    The energy terms are in k/2 form for bond, angle, dihedral and improper.
    The length is in unit of Angstrom, and angle is in unit of degree.
    The energies for bond, angle and vdW are in unit of kJ/mol/A^2, kJ/mol/rad^2 and kJ/mol.
    The two values for LJ interaction are sigma and epsilon.

    A force field term describing one topology element should only appear once.
    If there are duplicated terms, an Exception will be raised.

    Parameters
    ----------
    file : str

    Attributes
    ----------
    forcefield : ForceField
    '''

    def __init__(self, file):
        ff = ForceField()
        ff.lj_mixing_rule = ForceField.LJ_MIXING_GEOMETRIC
        ff.vdw_long_range = ForceField.VDW_LONGRANGE_CORRECT
        ff.vdw_cutoff = 1.2  # nm
        ff.scale_14_vdw = 0.5
        ff.scale_14_coulomb = 0.5

        self._parse(ff, file)
        self._setup(ff)

        self.forcefield = ff

    def _parse(self, ff, file):
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

            words = line.strip().split()
            if _section == 'atoms':
                self._parse_atom(ff, words)
            elif _section == 'bonds':
                self._parse_bond(ff, words)
            elif _section == 'angles':
                self._parse_angle(ff, words)
            elif _section == 'dihedrals':
                self._parse_dihedral(ff, words)
            elif _section == 'impropers':
                self._parse_improper(ff, words)
            elif _section == 'polarizations':
                self._parse_polarization(ff, words)
            else:
                raise Exception('Invalid header for fftool: %s' % line)

    def _setup(self, ff):
        # If this is a polarizable FF, add an AtomType and VdwTerm for Drude particles
        if ff.is_polarizable:
            type_drude = 'DP_'
            if type_drude not in ff.atom_types:
                dtype = AtomType(type_drude)
                vdw = LJ126Term(type_drude, type_drude, 0.0, 0.0)
                ff.add_term(dtype)
                ff.add_term(vdw)

        # If H* is defined in polar term, alpha of H will be merged into attached heavy atoms and the term will be removed
        if 'H*' in ff.polar_terms:
            hterm = ff.polar_terms.pop('H*')
            for pterm in ff.polar_terms.values():
                pterm.merge_alpha_H = hterm.alpha
            logger.info(f'H* found in polar term. '
                        f'Polarizability of H will be merged into attached heavy atoms')

    def _parse_atom(self, ff, words):
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
            ff.vdw_terms[lj.name] = lj
        else:
            raise Exception('Unsupported vdw potential form: %s' % (words[4]))

        if atype.name not in ff.atom_types.keys():
            ff.atom_types[atype.name] = atype
        else:
            raise Exception('Duplicated atom type: %s' % str(atype))

    def _parse_bond(self, ff, words):
        if words[2] not in ['harm', 'cons']:
            raise Exception('Unsupported bond function: %s' % (words[2]))
        term = HarmonicBondTerm(words[0], words[1], length=float(words[3]) / 10,
                                k=float(words[4]) * 100 / 2, fixed=(words[2] == 'cons'))
        if term.name in ff.bond_terms.keys():
            raise Exception('Duplicated bond term: %s' % str(term))
        ff.bond_terms[term.name] = term

    def _parse_angle(self, ff, words):
        if words[3] not in ['harm', 'cons']:
            raise Exception('Unsupported angle function: %s' % (words[3]))
        term = HarmonicAngleTerm(words[0], words[1], words[2], theta=float(words[4]) * DEG2RAD,
                                 k=float(words[5]) / 2, fixed=(words[3] == 'cons'))
        if term.name in ff.angle_terms.keys():
            raise Exception('Duplicated angle term: %s' % str(term))
        ff.angle_terms[term.name] = term

    def _parse_dihedral(self, ff, words):
        if words[4] not in ['opls']:
            raise Exception('Unsupported dihedral function: %s' % (words[4]))
        k1, k2, k3, k4 = list(map(lambda x: float(x) / 2, words[5:9]))
        term = OplsDihedralTerm(words[0], words[1], words[2], words[3], k1, k2, k3, k4)
        if term.name in ff.dihedral_terms.keys():
            raise Exception('Duplicated dihedral term: %s' % str(term))
        ff.dihedral_terms[term.name] = term

    def _parse_improper(self, ff, words):
        '''
        For improper i-j-k-l
        mstk always use CHARMM convention, i is the central atom
        fftool use OPLS convention, k is the central atom
        '''
        if words[4] not in ['opls']:
            raise Exception('Unsupported improper function: %s' % (words[4]))
        if {float(words[5]), float(words[7]), float(words[8])} != {0}:
            raise Exception('k1, k3 and k4 for opls improper should equal to 0')
        term = OplsImproperTerm(words[2], words[0], words[1], words[3], k=float(words[6]) / 2)
        if term.name in ff.improper_terms.keys():
            raise Exception('Duplicated improper term: %s' % str(term))
        ff.improper_terms[term.name] = term

    def _parse_polarization(slef, ff, words):
        name, mass, charge, k, alpha, thole = words
        term = DrudePolarTerm(name, float(alpha) / 1000, float(thole))
        term.mass = float(mass)
        term.k = float(k) / 2 * 100  # convert from kJ/mol/A^2 to kJ/mol/nm^2
        if term.name in ff.polar_terms.keys():
            raise Exception('Duplicated drude term: %s' % str(term))
        ff.polar_terms[term.name] = term


class Monomer:
    def __init__(self, name: str, charge: float, dipole: float, alpha: float = None):
        self.name = name
        self.charge = charge
        self.dipole = dipole
        self.alpha = alpha or -1  # set to -1 if not provided
        self.atoms: [str] = []  # list of atom types

    def __repr__(self):
        return '<Monomer: %s>' % self.name


class Dimer:
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


class PaduaLJScaler:
    '''
    PaduaLJScaler implements the empirical LJ scaling scheme of Goloviznina Kateryna.

    A scaling file should be providing,
    which provides information about fragments, dimers, atom types belongs to each fragments,
    charges and dipoles of each fragments and distances of dimers.

    Scaling factor can also be provided directly in DIMERS section.
    K_SAPT2 section is ignored.

    After the scaling file been loaded, all the LJ126Terms in a ForceField can be scaled by calling
    function :func:`scale`.

    Examples
    --------
    >>> scaler = PaduaLJScaler('ljscale.txt')
    >>> ff = ForceField.open('input.ff')
    >>> scaler.scale(ff)

    Parameters
    ----------
    file : str

    Attributes
    ----------
    scale_sigma : float
        Scaling factor for sigma, which equals to 1.0 by default.
        Other value can be provided by the input file by adding a line `SCALE_SIGMA value`

    '''

    def __init__(self, file):
        self.scale_sigma = 1.0
        self._monomers = []
        self._dimers = []
        self._parse(file)

    def _parse(self, file):
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
                self._monomers.append(Monomer(name, float(words[1]), float(words[2]), alpha))

            elif _section == 'DIMERS':
                m1 = next(m for m in self._monomers if m.name == words[0])
                m2 = next(m for m in self._monomers if m.name == words[1])
                scale = None
                if len(words) > 3:
                    try:
                        scale = float(words[3])
                    except:
                        pass
                self._dimers.append(Dimer(m1, m2, float(words[2]), scale))

            elif _section == 'ATOMS':
                monomer = next(m for m in self._monomers if m.name == words[0])
                monomer.atoms = words[1:]

    def scale(self, ff):
        '''
        Scale epsilon and sigma of all the LJ126Terms in a ForceField.

        Both the self-vdW terms and pairwise-vdW terms are scaled, as long as it is LJ126Term.

        Parameters
        ----------
        ff : ForceField

        Returns
        -------
        all_scaled : bool
            Whether or not all the vdW terms are been scaled.
            If it's False, it means fragment or dimer data not found for some LJ126Terms.
        '''
        _all_scaled = True
        for type1, type2 in itertools.combinations(ff.atom_types.values(), 2):
            # must scale pairwise vdW terms first
            # because they may be generated from self vdW terms by combination rule
            # the scaled pairwise terms must be added into the ff set
            try:
                vdw = ff.get_vdw_term(type1, type2)
            except FFTermNotFoundError:
                continue
            if type(vdw) is LJ126Term:
                if vdw.epsilon == 0:
                    continue
                if not self.scale_lj(vdw):
                    _all_scaled = False
                ff.add_term(vdw, replace=True)
        for vdw in ff.vdw_terms.values():
            if type(vdw) == LJ126Term:
                if vdw.epsilon == 0:
                    continue
                if not self.scale_lj(vdw):
                    _all_scaled = False

        return _all_scaled

    def scale_lj(self, term):
        '''
        Scale the epsilon and sigma of a LJ126Term.

        If scaling factor for epsilon is successfully predicted,
        then epsilon will be scaled by this factor, and sigma will be scaled by predefined scaling factor.
        If scaling factor for epsilon is not predicted (because of fragment or dimer not found),
        then both epsilon and sigma will be untouched.

        Parameters
        ----------
        term : LJ126Term

        Returns
        -------
        scaled : bool
        '''
        k_eps = self.predict_scale_epsilon(term.type1, term.type2)
        if k_eps is None:
            return False

        term.epsilon *= k_eps
        term.comments.append('eps*%.3f' % k_eps)
        k_sig = self.scale_sigma
        if k_sig != 1.0:
            term.sigma *= k_sig
            term.comments.append('sig*%.3f' % k_sig)

        return True

    def predict_scale_epsilon(self, at1, at2):
        '''
        Predict the scaling factor for epsilon of the LJ126Term between two atom types
        with empirical formulation.

        Parameters
        ----------
        at1 : str
            The first atom type
        at2 : str
            The second atom type

        Returns
        -------
        k : float or None
            float value means successful prediction.
            None means fragment or dimer data not found.
        '''
        m1 = m2 = None
        for monomer in self._monomers:
            if at1 in monomer.atoms:
                m1 = monomer
            if at2 in monomer.atoms:
                m2 = monomer
        if m1 is None or m2 is None:
            return None

        for dimer in self._dimers:
            if (m1 == dimer.monomer1 and m2 == dimer.monomer2) or (
                    m2 == dimer.monomer1 and m1 == dimer.monomer2):
                dimer = dimer
                break
        else:
            return None

        if dimer.scale_factor is None:
            dimer.predict_scale_epsilon()

        return dimer.scale_factor


ForceField.register_format('.ff', Padua)
