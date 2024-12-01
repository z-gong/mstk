import itertools
from enum import IntEnum
from mstk.simsys import System
from mstk.forcefield import *
from mstk.topology import *
from mstk.chem.constant import *

try:
    import openmm.openmm as mm
except ImportError:
    OPENMM_FOUND = False
else:
    OPENMM_FOUND = True


class ForceGroup(IntEnum):
    BOND = 1
    ANGLE = 2
    DIHEDRAL = 3
    IMPROPER = 4
    VDW = 5
    COULOMB = 6
    DRUDE = 7


class OpenMMExporter:
    '''
    OpenMMExporter export a :class:`System` to a OpenMM system.

    All of the force field terms supported by `mstk` can be exported to OpenMM.
    '''

    def __init__(self):
        pass

    @staticmethod
    def export(system, disable_inter_mol=False, **kwargs):
        '''
        Generate OpenMM system from a system

        Parameters
        ----------
        system : System
        disable_inter_mol : bool
            Disable inter-molecular interactions. This is useful for optimizing a batch of molecules

        Returns
        -------
        omm_system : openmm.openmm.System
        '''
        if not OPENMM_FOUND:
            raise ImportError('Can not import OpenMM')

        supported_terms = {LJ126Term, MieTerm, MorseVdwTerm,
                           HarmonicBondTerm, QuarticBondTerm, MorseBondTerm,
                           HarmonicAngleTerm, QuarticAngleTerm, HarmonicCosineAngleTerm, SDKAngleTerm, LinearAngleTerm,
                           OplsDihedralTerm, PeriodicDihedralTerm,
                           OplsImproperTerm, HarmonicImproperTerm,
                           DrudePolarTerm}
        unsupported = system.ff.energy_term_classes - supported_terms
        if unsupported != set():
            raise Exception('Unsupported FF terms: %s' % (', '.join(map(lambda x: x.__name__, unsupported))))

        if system.topology.virtual_site_classes - {TIP4PSite} != set():
            raise Exception('Virtual sites other than TIP4PSite haven\'t been implemented')

        omm_system = mm.System()
        if system.use_pbc:
            omm_system.setDefaultPeriodicBoxVectors(*system.topology.cell.vectors)
        for atom in system.topology.atoms:
            omm_system.addParticle(atom.mass)

        OpenMMExporter._setup_bonded(system, omm_system)
        if not disable_inter_mol:
            OpenMMExporter._setup_nonbonded(system, omm_system)
        else:
            OpenMMExporter._setup_nonbonded_as_bonded(system, omm_system)

        ### Remove COM motion ###################################################################
        logger.debug('Setting up COM motion remover...')
        omm_system.addForce(mm.CMMotionRemover(10))

        return omm_system

    @staticmethod
    def _setup_bonded(system, omm_system):
        top = system.topology
        ff = system.ff

        ### Set up bonds #######################################################################
        for bond_class in system.ff.bond_term_classes:
            if bond_class == HarmonicBondTerm:
                logger.debug('Setting up harmonic bonds...')
                bforce = mm.HarmonicBondForce()
                for bond in top.bonds:
                    if bond.is_drude:
                        # DrudeForce will handle the bond between Drude pair
                        continue
                    if bond in system.constrain_bonds:
                        continue
                    bterm = system.bond_terms[bond]
                    if type(bterm) != HarmonicBondTerm:
                        continue
                    bforce.addBond(bond.atom1.id, bond.atom2.id, bterm.length, bterm.k * 2)
            elif bond_class == QuarticBondTerm:
                logger.debug('Setting up quartic bonds...')
                bforce = mm.CustomBondForce('k2*(r-b0)^2+k3*(r-b0)^3+k4*(r-b0)^4')
                bforce.addPerBondParameter('b0')
                bforce.addPerBondParameter('k2')
                bforce.addPerBondParameter('k3')
                bforce.addPerBondParameter('k4')
                for bond in top.bonds:
                    if bond.is_drude:
                        # DrudeForce will handle the bond between Drude pair
                        continue
                    if bond in system.constrain_bonds:
                        continue
                    bterm = system.bond_terms[bond]
                    if type(bterm) != QuarticBondTerm:
                        continue
                    bforce.addBond(bond.atom1.id, bond.atom2.id, [bterm.length, bterm.k2, bterm.k3, bterm.k4])
            elif bond_class == MorseBondTerm:
                logger.debug('Setting up Morse bonds...')
                bforce = mm.CustomBondForce('depth*(1-exp(-alpha*(r-b0)))^2;'
                                            'alpha=(k/depth)^0.5')
                bforce.addPerBondParameter('b0')
                bforce.addPerBondParameter('k')
                bforce.addPerBondParameter('depth')
                for bond in top.bonds:
                    if bond.is_drude:
                        # DrudeForce will handle the bond between Drude pair
                        continue
                    if bond in system.constrain_bonds:
                        continue
                    bterm = system.bond_terms[bond]
                    if type(bterm) != MorseBondTerm:
                        continue
                    bforce.addBond(bond.atom1.id, bond.atom2.id, [bterm.length, bterm.k, bterm.depth])
            else:
                raise Exception('Bond terms other that HarmonicBondTerm, QuarticBondTerm and MorseBondTerm '
                                'haven\'t been implemented')
            bforce.setUsesPeriodicBoundaryConditions(system.use_pbc)
            bforce.setForceGroup(ForceGroup.BOND)
            bforce.setName('Bond')
            omm_system.addForce(bforce)

        ### Set up angles #######################################################################
        for angle_class in system.ff.angle_term_classes:
            if angle_class == HarmonicAngleTerm:
                logger.debug('Setting up harmonic angles...')
                aforce = mm.HarmonicAngleForce()
                for angle in top.angles:
                    if angle in system.constrain_angles:
                        continue
                    aterm = system.angle_terms[angle]
                    if type(aterm) == HarmonicAngleTerm:
                        aforce.addAngle(angle.atom1.id, angle.atom2.id, angle.atom3.id, aterm.theta, aterm.k * 2)
            elif angle_class == QuarticAngleTerm:
                logger.debug('Setting up quartic angles...')
                aforce = mm.CustomAngleForce('k2*(theta-theta0)^2+k3*(theta-theta0)^3+k4*(theta-theta0)^4')
                aforce.addPerAngleParameter('theta0')
                aforce.addPerAngleParameter('k2')
                aforce.addPerAngleParameter('k3')
                aforce.addPerAngleParameter('k4')
                for angle in top.angles:
                    if angle in system.constrain_angles:
                        continue
                    aterm = system.angle_terms[angle]
                    if type(aterm) == QuarticAngleTerm:
                        aforce.addAngle(angle.atom1.id, angle.atom2.id, angle.atom3.id,
                                        [aterm.theta, aterm.k2, aterm.k3, aterm.k4])
            elif angle_class == HarmonicCosineAngleTerm:
                logger.debug('Setting up harmonic cosine angles...')
                aforce = mm.CustomAngleForce('k/sin(theta0)^2*(cos(theta)-cos(theta0))^2')
                aforce.addPerAngleParameter('theta0')
                aforce.addPerAngleParameter('k')
                for angle in top.angles:
                    if angle in system.constrain_angles:
                        continue
                    aterm = system.angle_terms[angle]
                    if type(aterm) == HarmonicCosineAngleTerm:
                        aforce.addAngle(angle.atom1.id, angle.atom2.id, angle.atom3.id, [aterm.theta, aterm.k])
            elif angle_class == SDKAngleTerm:
                logger.debug('Setting up SDK angles...')
                aforce = mm.CustomCompoundBondForce(
                    3, 'k*(theta-theta0)^2+step(rmin-r)*LJ;'
                       'LJ=C*epsilon*((sigma/r)^n-(sigma/r)^m)+epsilon;'
                       'theta=angle(p1,p2,p3);'
                       'r=distance(p1,p3);'
                       'C=n/(n-m)*(n/m)^(m/(n-m));'
                       'rmin=(n/m)^(1/(n-m))*sigma')
                aforce.addPerBondParameter('theta0')
                aforce.addPerBondParameter('k')
                aforce.addPerBondParameter('epsilon')
                aforce.addPerBondParameter('sigma')
                aforce.addPerBondParameter('n')
                aforce.addPerBondParameter('m')
                for angle in top.angles:
                    if angle in system.constrain_angles:
                        continue
                    aterm = system.angle_terms[angle]
                    if type(aterm) != SDKAngleTerm:
                        continue
                    vdw = ff.get_vdw_term(ff.atom_types[angle.atom1.type], ff.atom_types[angle.atom3.type])
                    if type(vdw) != MieTerm:
                        raise Exception(f'Corresponding MieTerm for {aterm} not found in FF')
                    aforce.addBond([angle.atom1.id, angle.atom2.id, angle.atom3.id],
                                   [aterm.theta, aterm.k, vdw.epsilon, vdw.sigma, vdw.repulsion, vdw.attraction])
            elif angle_class == LinearAngleTerm:
                logger.debug('Setting up linear angles...')
                aforce = mm.CustomCompoundBondForce(
                    3, 'k_linear*r^2;'
                       'r=pointdistance(xref,yref,zref,x2,y2,z2);'
                       'xref=x1*a+x3*(1-a);'
                       'yref=y1*a+y3*(1-a);'
                       'zref=z1*a+z3*(1-a)')
                aforce.addPerBondParameter('a')
                aforce.addPerBondParameter('k_linear')
                for angle in top.angles:
                    if angle in system.constrain_angles:
                        continue
                    aterm = system.angle_terms[angle]
                    if type(aterm) != LinearAngleTerm:
                        continue
                    bond12, bond23 = angle.bonds
                    bterm12 = system.bond_terms[bond12]
                    bterm23 = system.bond_terms[bond23]
                    a, k_linear = aterm.calc_a_k_linear(bterm12.length, bterm23.length)
                    aforce.addBond([angle.atom1.id, angle.atom2.id, angle.atom3.id], [a, k_linear])
            else:
                raise Exception('Angle terms other that HarmonicAngleTerm, QuarticAngleTerm, HarmonicCosineAngleTerm, '
                                'SDKAngleTerm and LinearAngleTerm haven\'t been implemented')
            aforce.setUsesPeriodicBoundaryConditions(system.use_pbc)
            aforce.setForceGroup(ForceGroup.ANGLE)
            aforce.setName('Angle')
            omm_system.addForce(aforce)

        ### Set up constraints #################################################################
        logger.debug(f'Setting up {len(system.constrain_bonds)} bond constraints...')
        for bond in top.bonds:
            if bond in system.constrain_bonds:
                omm_system.addConstraint(bond.atom1.id, bond.atom2.id,
                                         system.constrain_bonds[bond])
        logger.debug(f'Setting up {len(system.constrain_angles)} angle constraints...')
        for angle in top.angles:
            if angle in system.constrain_angles:
                omm_system.addConstraint(angle.atom1.id, angle.atom3.id,
                                         system.constrain_angles[angle])

        ### Set up dihedrals ###################################################################
        for dihedral_class in system.ff.dihedral_term_classes:
            if dihedral_class in (OplsDihedralTerm, PeriodicDihedralTerm):
                logger.debug('Setting up periodic dihedrals...')
                dforce = mm.PeriodicTorsionForce()
                for dihedral in top.dihedrals:
                    dterm = system.dihedral_terms[dihedral]
                    ia1, ia2, ia3, ia4 = dihedral.atom1.id, dihedral.atom2.id, dihedral.atom3.id, dihedral.atom4.id
                    if type(dterm) is OplsDihedralTerm:
                        dterm = dterm.to_periodic_term()
                    if type(dterm) == PeriodicDihedralTerm:
                        for phase in dterm.phases:
                            if phase.k != 0:
                                dforce.addTorsion(ia1, ia2, ia3, ia4, phase.n, phase.phi, phase.k)
            else:
                raise Exception('Dihedral terms other that OplsDihedralTerm and PeriodicDihedralTerm '
                                'haven\'t been implemented')
            dforce.setUsesPeriodicBoundaryConditions(system.use_pbc)
            dforce.setForceGroup(ForceGroup.DIHEDRAL)
            dforce.setName('Dihedral')
            omm_system.addForce(dforce)

        ### Set up impropers ####################################################################
        for improper_class in system.ff.improper_term_classes:
            if improper_class == OplsImproperTerm:
                logger.debug('Setting up periodic impropers...')
                iforce = mm.CustomTorsionForce('k*(1-cos(2*theta))')
                iforce.addPerTorsionParameter('k')
                for improper in top.impropers:
                    iterm = system.improper_terms[improper]
                    if type(iterm) == OplsImproperTerm:
                        # in OPLS convention, the third atom is the central atom
                        iforce.addTorsion(improper.atom2.id, improper.atom3.id,
                                          improper.atom1.id, improper.atom4.id, [iterm.k])
            elif improper_class == HarmonicImproperTerm:
                logger.debug('Setting up harmonic impropers...')
                iforce = mm.CustomTorsionForce(f'k*min(dtheta,2*{PI}-dtheta)^2;'
                                               f'dtheta=abs(theta-phi0)')
                iforce.addPerTorsionParameter('phi0')
                iforce.addPerTorsionParameter('k')
                for improper in top.impropers:
                    iterm = system.improper_terms[improper]
                    if type(iterm) == HarmonicImproperTerm:
                        iforce.addTorsion(improper.atom1.id, improper.atom2.id,
                                          improper.atom3.id, improper.atom4.id,
                                          [iterm.phi, iterm.k])
            else:
                raise Exception('Improper terms other that PeriodicImproperTerm and '
                                'HarmonicImproperTerm haven\'t been implemented')
            iforce.setUsesPeriodicBoundaryConditions(system.use_pbc)
            iforce.setForceGroup(ForceGroup.IMPROPER)
            iforce.setName('Improper')
            omm_system.addForce(iforce)

    @staticmethod
    def _setup_nonbonded(system, omm_system):
        top = system.topology
        ff = system.ff

        ### Set up non-bonded interactions #########################################################
        # NonbonedForce is not flexible enough. Use it only for Coulomb interactions (including 1-4 Coulomb exceptions)
        # CustomNonbondedForce handles vdW interactions
        # CustomBondForce handles 1-4 vdW exceptions
        cutoff = ff.vdw_cutoff
        logger.debug('Setting up Coulomb interactions...')
        nbforce = mm.NonbondedForce()
        if system.use_pbc:
            nbforce.setNonbondedMethod(mm.NonbondedForce.PME)
            nbforce.setEwaldErrorTolerance(5E-4)
            nbforce.setCutoffDistance(cutoff)
            # dispersion will be handled by CustomNonbondedForce
            nbforce.setUseDispersionCorrection(False)
            nbforce.setExceptionsUsePeriodicBoundaryConditions(True)
        else:
            nbforce.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        nbforce.setForceGroup(ForceGroup.COULOMB)
        nbforce.setName('Coulomb')
        for atom in top.atoms:
            nbforce.addParticle(atom.charge, 1.0, 0.0)
        if system.charged:
            omm_system.addForce(nbforce)

        ### Set up vdW interactions #########################################################
        atom_types = list(ff.atom_types.values())
        type_names = list(ff.atom_types.keys())
        n_type = len(atom_types)
        for vdw_class in system.ff.vdw_term_classes:
            if vdw_class == LJ126Term:
                logger.debug('Setting up LJ-12-6 vdW interactions...')
                if system.use_pbc and ff.vdw_long_range == ForceField.VDW_LONGRANGE_SHIFT:
                    invRc6 = 1 / cutoff ** 6
                    cforce = mm.CustomNonbondedForce(
                        f'A(type1,type2)*(invR6*invR6-{invRc6 * invRc6})-'
                        f'B(type1,type2)*(invR6-{invRc6});'
                        f'invR6=1/r^6')
                else:
                    cforce = mm.CustomNonbondedForce(
                        'A(type1,type2)*invR6*invR6-B(type1,type2)*invR6;'
                        'invR6=1/r^6')
                cforce.addPerParticleParameter('type')
                A_list = [0.0] * n_type * n_type
                B_list = [0.0] * n_type * n_type
                for i, atype1 in enumerate(atom_types):
                    for j, atype2 in enumerate(atom_types):
                        vdw = ff.get_vdw_term(atype1, atype2)
                        if type(vdw) == LJ126Term:
                            A = 4 * vdw.epsilon * vdw.sigma ** 12
                            B = 4 * vdw.epsilon * vdw.sigma ** 6
                        else:
                            A = B = 0
                        A_list[i + n_type * j] = A
                        B_list[i + n_type * j] = B
                cforce.addTabulatedFunction('A', mm.Discrete2DFunction(n_type, n_type, A_list))
                cforce.addTabulatedFunction('B', mm.Discrete2DFunction(n_type, n_type, B_list))

                for atom in top.atoms:
                    id_type = type_names.index(atom.type)
                    cforce.addParticle([id_type])

            elif vdw_class == MorseVdwTerm:
                logger.debug('Setting up Morse vdW interactions...')
                if system.use_pbc and ff.vdw_long_range == ForceField.VDW_LONGRANGE_SHIFT:
                    cforce = mm.CustomNonbondedForce(
                        f'DEPTH(type1,type2)*(1-exp(-ALPHA(type1,type2)*(r-R0(type1,type2))))^2-'
                        f'DEPTH(type1,type2)*(1-exp(-ALPHA(type1,type2)*({cutoff}-R0(type1,type2))))^2'
                    )
                else:
                    cforce = mm.CustomNonbondedForce(
                        f'DEPTH(type1,type2)*((1-exp(-ALPHA(type1,type2)*(r-R0(type1,type2))))^2-1)'
                    )
                cforce.addPerParticleParameter('type')
                DEPTH_list = [0.0] * n_type * n_type
                ALPHA_list = [0.0] * n_type * n_type
                R0_list = [0.0] * n_type * n_type
                for i, atype1 in enumerate(atom_types):
                    for j, atype2 in enumerate(atom_types):
                        vdw = ff.get_vdw_term(atype1, atype2)
                        if type(vdw) == MorseVdwTerm:
                            DEPTH = vdw.depth
                            ALPHA = vdw.alpha
                            R0 = vdw.r0
                        else:
                            DEPTH = ALPHA = R0 = 0
                        DEPTH_list[i + n_type * j] = DEPTH
                        ALPHA_list[i + n_type * j] = ALPHA
                        R0_list[i + n_type * j] = R0
                cforce.addTabulatedFunction('DEPTH', mm.Discrete2DFunction(n_type, n_type, DEPTH_list))
                cforce.addTabulatedFunction('ALPHA', mm.Discrete2DFunction(n_type, n_type, ALPHA_list))
                cforce.addTabulatedFunction('R0', mm.Discrete2DFunction(n_type, n_type, R0_list))

                for atom in top.atoms:
                    id_type = type_names.index(atom.type)
                    cforce.addParticle([id_type])

            elif vdw_class == MieTerm:
                logger.debug('Setting up Mie vdW interactions...')
                if system.use_pbc and ff.vdw_long_range == ForceField.VDW_LONGRANGE_SHIFT:
                    cforce = mm.CustomNonbondedForce('A(type1,type2)/r^REP(type1,type2)-'
                                                     'B(type1,type2)/r^ATT(type1,type2)-'
                                                     'SHIFT(type1,type2)')
                else:
                    cforce = mm.CustomNonbondedForce('A(type1,type2)/r^REP(type1,type2)-'
                                                     'B(type1,type2)/r^ATT(type1,type2)')
                cforce.addPerParticleParameter('type')
                A_list = [0.0] * n_type * n_type
                B_list = [0.0] * n_type * n_type
                REP_list = [0.0] * n_type * n_type
                ATT_list = [0.0] * n_type * n_type
                SHIFT_list = [0.0] * n_type * n_type
                for i, atype1 in enumerate(atom_types):
                    for j, atype2 in enumerate(atom_types):
                        vdw = ff.get_vdw_term(atype1, atype2)
                        if type(vdw) == MieTerm:
                            A = vdw.factor_energy() * vdw.epsilon * vdw.sigma ** vdw.repulsion
                            B = vdw.factor_energy() * vdw.epsilon * vdw.sigma ** vdw.attraction
                            REP = vdw.repulsion
                            ATT = vdw.attraction
                            SHIFT = A / cutoff ** REP - B / cutoff ** ATT
                        else:
                            A = B = REP = ATT = SHIFT = 0
                        A_list[i + n_type * j] = A
                        B_list[i + n_type * j] = B
                        REP_list[i + n_type * j] = REP
                        ATT_list[i + n_type * j] = ATT
                        SHIFT_list[i + n_type * j] = SHIFT
                cforce.addTabulatedFunction('A', mm.Discrete2DFunction(n_type, n_type, A_list))
                cforce.addTabulatedFunction('B', mm.Discrete2DFunction(n_type, n_type, B_list))
                cforce.addTabulatedFunction('REP', mm.Discrete2DFunction(n_type, n_type, REP_list))
                cforce.addTabulatedFunction('ATT', mm.Discrete2DFunction(n_type, n_type, ATT_list))
                if system.use_pbc and ff.vdw_long_range == ForceField.VDW_LONGRANGE_SHIFT:
                    cforce.addTabulatedFunction('SHIFT', mm.Discrete2DFunction(n_type, n_type, SHIFT_list))

                for atom in top.atoms:
                    id_type = type_names.index(atom.type)
                    cforce.addParticle([id_type])

            else:
                raise Exception('vdW terms other than LJ126Term, MieTerm and MorseVdwTerm '
                                'haven\'t been implemented')
            if system.use_pbc:
                cforce.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
                cforce.setCutoffDistance(cutoff)
                if ff.vdw_long_range == ForceField.VDW_LONGRANGE_CORRECT:
                    cforce.setUseLongRangeCorrection(True)
            else:
                cforce.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
            cforce.setForceGroup(ForceGroup.VDW)
            cforce.setName('vdW')
            omm_system.addForce(cforce)

        ### Set up 1-2, 1-3 and 1-4 exceptions ##################################################
        logger.debug('Setting up 1-2, 1-3 and 1-4 exceptions...')
        custom_nb_forces = [f for f in omm_system.getForces() if type(f) == mm.CustomNonbondedForce]
        pair12, pair13, pair14 = top.get_12_13_14_pairs()
        for atom1, atom2 in pair12 + pair13:
            nbforce.addException(atom1.id, atom2.id, 0.0, 1.0, 0.0)
            for f in custom_nb_forces:
                f.addExclusion(atom1.id, atom2.id)
        # As long as 1-4 LJ OR Coulomb need to be scaled, then this pair should be excluded from ALL non-bonded forces.
        # This is required by OpenMM's internal implementation.
        # Even though NonbondedForce can handle 1-4 vdW, we use it only for 1-4 Coulomb.
        # And use CustomBondForce to handle 1-4 vdW, which makes it more clear for energy decomposition.
        if ff.scale_14_vdw != 1 or ff.scale_14_coulomb != 1:
            pair14_forces = {}  # {VdwTerm: mm.NbForce}
            for atom1, atom2 in pair14:
                charge_prod = atom1.charge * atom2.charge * ff.scale_14_coulomb
                nbforce.addException(atom1.id, atom2.id, charge_prod, 1.0, 0.0)
                for f in custom_nb_forces:
                    f.addExclusion(atom1.id, atom2.id)
                if ff.scale_14_vdw == 0:
                    continue
                vdw = ff.get_vdw_term(ff.atom_types[atom1.type], ff.atom_types[atom2.type])
                # We generalize LJ126Term and MieTerm because of minimal computational cost for 1-4 vdW
                if type(vdw) in (LJ126Term, MieTerm):
                    cbforce = pair14_forces.get(MieTerm)
                    if cbforce is None:
                        cbforce = mm.CustomBondForce('C*epsilon*((sigma/r)^n-(sigma/r)^m);'
                                                     'C=n/(n-m)*(n/m)^(m/(n-m))')
                        cbforce.addPerBondParameter('epsilon')
                        cbforce.addPerBondParameter('sigma')
                        cbforce.addPerBondParameter('n')
                        cbforce.addPerBondParameter('m')
                        cbforce.setUsesPeriodicBoundaryConditions(system.use_pbc)
                        cbforce.setForceGroup(ForceGroup.VDW)
                        cbforce.setName('vdW')
                        omm_system.addForce(cbforce)
                        pair14_forces[MieTerm] = cbforce
                    epsilon = vdw.epsilon * ff.scale_14_vdw
                    if type(vdw) == LJ126Term:
                        cbforce.addBond(atom1.id, atom2.id, [epsilon, vdw.sigma, 12, 6])
                    elif type(vdw) == MieTerm:
                        cbforce.addBond(atom1.id, atom2.id,
                                        [epsilon, vdw.sigma, vdw.repulsion, vdw.attraction])
                elif type(vdw) is MorseVdwTerm:
                    cbforce = pair14_forces.get(MorseVdwTerm)
                    if cbforce is None:
                        cbforce = mm.CustomBondForce('DEPTH*((1-exp(-ALPHA*(r-R0)))^2-1)')
                        cbforce.addPerBondParameter('DEPTH')
                        cbforce.addPerBondParameter('ALPHA')
                        cbforce.addPerBondParameter('R0')
                        cbforce.setUsesPeriodicBoundaryConditions(system.use_pbc)
                        cbforce.setForceGroup(ForceGroup.VDW)
                        cbforce.setName('vdW')
                        omm_system.addForce(cbforce)
                        pair14_forces[MorseVdwTerm] = cbforce
                    depth = vdw.depth * ff.scale_14_vdw
                    cbforce.addBond(atom1.id, atom2.id, [depth, vdw.alpha, vdw.r0])
                else:
                    raise Exception('1-4 scaling for vdW terms other than LJ126Term, MieTerm and MorseVdwTerm '
                                    'haven\'t been implemented')

        ### Set up Drude particles ##############################################################
        for polar_class in system.ff.polar_term_classes:
            if polar_class == DrudePolarTerm:
                logger.debug('Setting up Drude polarizations...')
                pforce = mm.DrudeForce()
                pforce.setForceGroup(ForceGroup.DRUDE)
                pforce.setName('Drude')
                omm_system.addForce(pforce)
                parent_idx_thole = {}  # {parent: (index in DrudeForce, thole)} for addScreenPair
                for parent, drude in system.drude_pairs.items():
                    pterm = system.polar_terms[parent]
                    n_H = len([atom for atom in parent.bond_partners if atom.symbol == 'H'])
                    alpha = pterm.alpha + n_H * pterm.merge_alpha_H
                    idx = pforce.addParticle(drude.id, parent.id, -1, -1, -1,
                                             drude.charge, alpha, 0, 0)
                    parent_idx_thole[parent] = (idx, pterm.thole)

                # exclude the non-boned interactions between Drude and parent
                # and those concerning Drude particles in 1-2 and 1-3 pairs
                # pairs formed by real atoms have already been handled above
                # also apply thole screening between 1-2 and 1-3 Drude dipole pairs
                drude_exclusions = list(system.drude_pairs.items())
                for atom1, atom2 in pair12 + pair13:
                    drude1 = system.drude_pairs.get(atom1)
                    drude2 = system.drude_pairs.get(atom2)
                    if drude1 is not None:
                        drude_exclusions.append((drude1, atom2))
                    if drude2 is not None:
                        drude_exclusions.append((atom1, drude2))
                    if drude1 is not None and drude2 is not None:
                        drude_exclusions.append((drude1, drude2))
                        idx1, thole1 = parent_idx_thole[atom1]
                        idx2, thole2 = parent_idx_thole[atom2]
                        pforce.addScreenedPair(idx1, idx2, (thole1 + thole2) / 2)
                for a1, a2 in drude_exclusions:
                    nbforce.addException(a1.id, a2.id, 0, 1.0, 0)
                    for f in custom_nb_forces:
                        f.addExclusion(a1.id, a2.id)

                # scale the non-boned interactions concerning Drude particles in 1-4 pairs
                # pairs formed by real atoms have already been handled above
                drude_exceptions14 = []
                for atom1, atom2 in pair14:
                    drude1 = system.drude_pairs.get(atom1)
                    drude2 = system.drude_pairs.get(atom2)
                    if drude1 is not None:
                        drude_exceptions14.append((drude1, atom2))
                    if drude2 is not None:
                        drude_exceptions14.append((atom1, drude2))
                    if drude1 is not None and drude2 is not None:
                        drude_exceptions14.append((drude1, drude2))
                for a1, a2 in drude_exceptions14:
                    charge_prod = a1.charge * a2.charge * ff.scale_14_coulomb
                    nbforce.addException(a1.id, a2.id, charge_prod, 1.0, 0.0)
                    for f in custom_nb_forces:
                        f.addExclusion(a1.id, a2.id)
            else:
                raise Exception('Polar terms other that DrudePolarTerm haven\'t been implemented')

        ### Set up virtual sites ################################################################
        if top.has_virtual_site:
            logger.debug('Setting up virtual sites...')
            for atom in top.atoms:
                vsite = atom.virtual_site
                if type(vsite) == TIP4PSite:
                    O, H1, H2 = vsite.parents
                    coeffs = system.get_TIP4P_linear_coeffs(atom)
                    omm_vsite = mm.ThreeParticleAverageSite(O.id, H1.id, H2.id, *coeffs)
                    omm_system.setVirtualSite(atom.id, omm_vsite)
                elif vsite is not None:
                    raise Exception('Virtual sites other than TIP4PSite haven\'t been implemented')

            # exclude the non-boned interactions between virtual sites and parents
            # and particles (atoms, drude particles, virtual sites) in 1-2 and 1-3 pairs
            # TODO Assume no more than one virtual site is attached to each atom
            vsite_exclusions = list(system.vsite_pairs.items())
            for atom, vsite in system.vsite_pairs.items():
                drude = system.drude_pairs.get(atom)
                if drude is not None:
                    vsite_exclusions.append((vsite, drude))
            for atom1, atom2 in pair12 + pair13:
                vsite1 = system.vsite_pairs.get(atom1)
                vsite2 = system.vsite_pairs.get(atom2)
                drude1 = system.drude_pairs.get(atom1)
                drude2 = system.drude_pairs.get(atom2)
                if vsite1 is not None:
                    vsite_exclusions.append((vsite1, atom2))
                    if drude2 is not None:
                        vsite_exclusions.append((vsite1, drude2))
                if vsite2 is not None:
                    vsite_exclusions.append((vsite2, atom1))
                    if drude1 is not None:
                        vsite_exclusions.append((vsite2, drude1))
                if None not in [vsite1, vsite2]:
                    vsite_exclusions.append((vsite1, vsite2))
            for a1, a2 in vsite_exclusions:
                nbforce.addException(a1.id, a2.id, 0, 1.0, 0)
                for f in custom_nb_forces:
                    f.addExclusion(a1.id, a2.id)

            # scale the non-boned interactions between virtual sites and particles in 1-4 pairs
            # TODO Assume no 1-4 LJ interactions on virtual sites
            vsite_exceptions14 = []
            for atom1, atom2 in pair14:
                vsite1 = system.vsite_pairs.get(atom1)
                vsite2 = system.vsite_pairs.get(atom2)
                drude1 = system.drude_pairs.get(atom1)
                drude2 = system.drude_pairs.get(atom2)
                if vsite1 is not None:
                    vsite_exceptions14.append((vsite1, atom2))
                    if drude2 is not None:
                        vsite_exceptions14.append((vsite1, drude2))
                if vsite2 is not None:
                    vsite_exceptions14.append((vsite2, atom1))
                    if drude1 is not None:
                        vsite_exceptions14.append((vsite2, drude1))
                if None not in [vsite1, vsite2]:
                    vsite_exceptions14.append((vsite1, vsite2))
            for a1, a2 in vsite_exceptions14:
                charge_prod = a1.charge * a2.charge * ff.scale_14_coulomb
                nbforce.addException(a1.id, a2.id, charge_prod, 1.0, 0.0)
                for f in custom_nb_forces:
                    f.addExclusion(a1.id, a2.id)

    @staticmethod
    def _setup_nonbonded_as_bonded(system, omm_system):
        if system.ff.polar_term_classes:
            raise Exception('Polar terms not supported')

        if system.topology.virtual_site_classes:
            raise Exception('Virtual sites not supported')

        if system.use_pbc:
            raise Exception('PBC not supported')

        top = system.topology
        ff = system.ff

        logger.debug('Setting up Coulomb interactions...')
        qqforce = mm.CustomBondForce(f'{ONE_4PI_EPS0}*qq/r')
        qqforce.addPerBondParameter('qq')
        qqforce.setUsesPeriodicBoundaryConditions(system.use_pbc)
        qqforce.setForceGroup(ForceGroup.COULOMB)
        qqforce.setName('Coulomb')
        omm_system.addForce(qqforce)

        logger.debug('Setting up vdW interactions...')
        ljforce = mm.CustomBondForce('C*epsilon*((sigma/r)^n-(sigma/r)^m);'
                                     'C=n/(n-m)*(n/m)^(m/(n-m))')
        ljforce.addPerBondParameter('epsilon')
        ljforce.addPerBondParameter('sigma')
        ljforce.addPerBondParameter('n')
        ljforce.addPerBondParameter('m')
        ljforce.setUsesPeriodicBoundaryConditions(system.use_pbc)
        ljforce.setForceGroup(ForceGroup.VDW)
        ljforce.setName('vdW')
        omm_system.addForce(ljforce)

        for mol in top.molecules:
            pairs12, pairs13, pairs14 = mol.get_12_13_14_pairs()
            # each pair tuple (Atom, Atom) in each pairs list are sorted by the Atom.id attributes of two atoms
            pairs1n = list(set(itertools.combinations(mol.atoms, 2)) - set(pairs12).union(pairs13).union(pairs14))
            for i, (atom1, atom2) in enumerate(pairs14 + pairs1n):
                vdw = ff.get_vdw_term(ff.atom_types[atom1.type], ff.atom_types[atom2.type])
                epsilon = vdw.epsilon
                qq = atom1.charge * atom2.charge
                if i < len(pairs14):
                    epsilon *= ff.scale_14_vdw
                    qq *= ff.scale_14_coulomb
                if type(vdw) == LJ126Term:
                    ljforce.addBond(atom1.id, atom2.id, [epsilon, vdw.sigma, 12, 6])
                elif type(vdw) == MieTerm:
                    ljforce.addBond(atom1.id, atom2.id, [epsilon, vdw.sigma, vdw.repulsion, vdw.attraction])
                else:
                    raise Exception(f'Cannot setup {vdw} as a bonded term')
                qqforce.addBond(atom1.id, atom2.id, [qq])
