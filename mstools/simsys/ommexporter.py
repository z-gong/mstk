from enum import IntEnum
from .system import System
from ..forcefield import *
from ..topology import *
from .. import logger

class ForceGroup(IntEnum):
    BOND = 1
    ANGLE = 2
    DIHEDRAL = 3
    IMPROPER = 4
    VDW = 5
    COULOMB = 6
    DRUDE = 7

class OpenMMExporter():
    '''
    OpenMMExporter export a :class:`System` to a OpenMM system
    '''
    def __init__(self):
        pass

    @staticmethod
    def export(system):
        '''
        Generate OpenMM system from a system

        Parameters
        ----------
        system : System

        Returns
        -------
        omm_system : simtk.openmm.System
        '''
        try:
            import simtk.openmm as mm
        except ImportError:
            raise ImportError('Can not import OpenMM')

        supported_terms = {LJ126Term, MieTerm,
                           HarmonicBondTerm,
                           HarmonicAngleTerm, SDKAngleTerm,
                           PeriodicDihedralTerm,
                           OplsImproperTerm, HarmonicImproperTerm,
                           DrudeTerm}
        unsupported = system.ff_classes - supported_terms
        if unsupported != set():
            raise Exception('Unsupported FF terms: %s'
                            % (', '.join(map(lambda x: x.__name__, unsupported))))

        omm_system = mm.System()
        omm_system.setDefaultPeriodicBoxVectors(*system._topology.cell.vectors)
        for atom in system._topology.atoms:
            omm_system.addParticle(atom.mass)

        ### Set up bonds #######################################################################
        for bond_class in system.bond_classes:
            if bond_class == HarmonicBondTerm:
                logger.info('Setting up harmonic bonds...')
                bforce = mm.HarmonicBondForce()
                for bond in system._topology.bonds:
                    if bond.is_drude:
                        continue
                    bterm = system.bond_terms[id(bond)]
                    if type(bterm) != HarmonicBondTerm:
                        continue
                    bforce.addBond(bond.atom1.id, bond.atom2.id, bterm.length, bterm.k * 2)
            else:
                raise Exception('Bond terms other that HarmonicBondTerm '
                                'haven\'t been implemented')
            bforce.setUsesPeriodicBoundaryConditions(True)
            bforce.setForceGroup(ForceGroup.BOND)
            omm_system.addForce(bforce)

        ### Set up angles #######################################################################
        for angle_class in system.angle_classes:
            if angle_class == HarmonicAngleTerm:
                logger.info('Setting up harmonic angles...')
                aforce = mm.HarmonicAngleForce()
                for angle in system._topology.angles:
                    aterm = system.angle_terms[id(angle)]
                    if type(aterm) == HarmonicAngleTerm:
                        aforce.addAngle(angle.atom1.id, angle.atom2.id, angle.atom3.id,
                                        aterm.theta * PI / 180, aterm.k * 2)
            elif angle_class == SDKAngleTerm:
                logger.info('Setting up SDK angles...')
                aforce = mm.CustomCompoundBondForce(
                    3, 'k*(theta-theta0)^2+step(rmin-r)*LJ96;'
                       'LJ96=6.75*epsilon*((sigma/r)^9-(sigma/r)^6)+epsilon;'
                       'theta=angle(p1,p2,p3);'
                       'r=distance(p1,p3);'
                       'rmin=1.144714*sigma')
                aforce.addPerBondParameter('theta0')
                aforce.addPerBondParameter('k')
                aforce.addPerBondParameter('epsilon')
                aforce.addPerBondParameter('sigma')
                for angle in system._topology.angles:
                    aterm = system.angle_terms[id(angle)]
                    if type(aterm) != SDKAngleTerm:
                        continue
                    vdw = system._ff.get_vdw_term(angle.atom1.type, angle.atom2.type)
                    if type(vdw) != MieTerm or vdw.repulsion != 9 or vdw.attraction != 6:
                        raise Exception(f'Corresponding 9-6 MieTerm for {aterm} not found in FF')
                    aforce.addBond([angle.atom1.id, angle.atom2.id, angle.atom3.id],
                                   [aterm.theta * PI / 180, aterm.k, vdw.epsilon, vdw.sigma])
            else:
                raise Exception('Angle terms other that HarmonicAngleTerm and SDKAngleTerm '
                                'haven\'t been implemented')
            aforce.setUsesPeriodicBoundaryConditions(True)
            aforce.setForceGroup(ForceGroup.ANGLE)
            omm_system.addForce(aforce)

        ### Set up constraints #################################################################
        logger.info(f'Setting up {len(system.constrain_bonds)} bond constraints...')
        for bond in system._topology.bonds:
            if id(bond) in system.constrain_bonds:
                omm_system.addConstraint(bond.atom1.id, bond.atom2.id,
                                         system.constrain_bonds[id(bond)])
        logger.info(f'Setting up {len(system.constrain_angles)} angle constraints...')
        for angle in system._topology.angles:
            if id(angle) in system.constrain_angles:
                omm_system.addConstraint(angle.atom1.id, angle.atom3.id,
                                         system.constrain_angles[id(angle)])

        ### Set up dihedrals ###################################################################
        for dihedral_class in system.dihedral_classes:
            if dihedral_class == PeriodicDihedralTerm:
                logger.info('Setting up periodic dihedrals...')
                dforce = mm.PeriodicTorsionForce()
                for dihedral in system._topology.dihedrals:
                    dterm = system.dihedral_terms[id(dihedral)]
                    ia1, ia2, ia3, ia4 = dihedral.atom1.id, dihedral.atom2.id, dihedral.atom3.id, dihedral.atom4.id
                    if type(dterm) == PeriodicDihedralTerm:
                        for par in dterm.parameters:
                            dforce.addTorsion(ia1, ia2, ia3, ia4, par.n, par.phi * PI / 180, par.k)
                    else:
                        continue
            else:
                raise Exception('Dihedral terms other that PeriodicDihedralTerm '
                                'haven\'t been implemented')
            dforce.setUsesPeriodicBoundaryConditions(True)
            dforce.setForceGroup(ForceGroup.DIHEDRAL)
            omm_system.addForce(dforce)

        ### Set up impropers ####################################################################
        for improper_class in system.improper_classes:
            if improper_class == OplsImproperTerm:
                logger.info('Setting up periodic impropers...')
                iforce = mm.CustomTorsionForce('k*(1-cos(2*theta))')
                iforce.addPerTorsionParameter('k')
                for improper in system._topology.impropers:
                    iterm = system.improper_terms[id(improper)]
                    if type(iterm) == OplsImproperTerm:
                        # in OPLS convention, the third atom is the central atom
                        iforce.addTorsion(improper.atom2.id, improper.atom3.id,
                                          improper.atom1.id, improper.atom4.id, [iterm.k])
            elif improper_class == HarmonicImproperTerm:
                logger.info('Setting up harmonic impropers...')
                iforce = mm.CustomTorsionForce(f'k*min(dtheta,2*pi-dtheta)^2;'
                                               f'dtheta=abs(theta-phi0);'
                                               f'pi={PI}')
                iforce.addPerTorsionParameter('phi0')
                iforce.addPerTorsionParameter('k')
                for improper in system._topology.impropers:
                    iterm = system.improper_terms[id(improper)]
                    if type(iterm) == HarmonicImproperTerm:
                        iforce.addTorsion(improper.atom1.id, improper.atom2.id,
                                          improper.atom3.id, improper.atom4.id,
                                          [iterm.phi * PI / 180, iterm.k])
            else:
                raise Exception('Improper terms other that PeriodicImproperTerm and '
                                'HarmonicImproperTerm haven\'t been implemented')
            iforce.setUsesPeriodicBoundaryConditions(True)
            iforce.setForceGroup(ForceGroup.IMPROPER)
            omm_system.addForce(iforce)

        ### Set up non-bonded interactions #########################################################
        # NonbonedForce is not flexible enough. Use it only for Coulomb interactions
        # CustomNonbondedForce handles vdW interactions
        cutoff = system._ff.vdw_cutoff
        logger.info('Setting up Coulomb interactions...')
        nbforce = mm.NonbondedForce()
        nbforce.setNonbondedMethod(mm.NonbondedForce.PME)
        nbforce.setEwaldErrorTolerance(5E-4)
        nbforce.setCutoffDistance(cutoff)
        nbforce.setUseDispersionCorrection(False)
        nbforce.setForceGroup(ForceGroup.COULOMB)
        omm_system.addForce(nbforce)
        for atom in system._topology.atoms:
            nbforce.addParticle(atom.charge, 1.0, 0.0)

        ### Set up vdW interactions #########################################################
        atom_types = list(system._ff.atom_types.values())
        type_names = list(system._ff.atom_types.keys())
        n_type = len(atom_types)
        for vdw_class in system.vdw_classes:
            if vdw_class == LJ126Term:
                logger.info('Setting up LJ-12-6 vdW interactions...')
                if system._ff.vdw_long_range == ForceField.VDW_LONGRANGE_SHIFT:
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
                for i, type1 in enumerate(atom_types):
                    for j, type2 in enumerate(atom_types):
                        vdw = system._ff.get_vdw_term(type1, type2)
                        if type(vdw) == LJ126Term:
                            A = 4 * vdw.epsilon * vdw.sigma ** 12
                            B = 4 * vdw.epsilon * vdw.sigma ** 6
                        else:
                            A = B = 0
                        A_list[i + n_type * j] = A
                        B_list[i + n_type * j] = B
                cforce.addTabulatedFunction('A', mm.Discrete2DFunction(n_type, n_type, A_list))
                cforce.addTabulatedFunction('B', mm.Discrete2DFunction(n_type, n_type, B_list))

                for atom in system._topology.atoms:
                    id_type = type_names.index(atom.type)
                    cforce.addParticle([id_type])

            elif vdw_class == MieTerm:
                logger.info('Setting up Mie vdW interactions...')
                if system._ff.vdw_long_range == ForceField.VDW_LONGRANGE_SHIFT:
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
                for i, type1 in enumerate(atom_types):
                    for j, type2 in enumerate(atom_types):
                        vdw = system._ff.get_vdw_term(type1, type2)
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
                if system._ff.vdw_long_range == ForceField.VDW_LONGRANGE_SHIFT:
                    cforce.addTabulatedFunction('SHIFT',
                                                mm.Discrete2DFunction(n_type, n_type, SHIFT_list))

                for atom in system._topology.atoms:
                    id_type = type_names.index(atom.type)
                    cforce.addParticle([id_type])

            else:
                raise Exception('vdW terms other than LJ126Term and MieTerm '
                                'haven\'t been implemented')
            cforce.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            cforce.setCutoffDistance(cutoff)
            if system._ff.vdw_long_range == ForceField.VDW_LONGRANGE_CORRECT:
                cforce.setUseLongRangeCorrection(True)
            else:
                cforce.setUseLongRangeCorrection(False)
            cforce.setForceGroup(ForceGroup.VDW)
            omm_system.addForce(cforce)

        ### Set up 1-2, 1-3 and 1-4 exceptions ##################################################
        logger.info('Setting up 1-2, 1-3 and 1-4 vdW interactions...')
        custom_nb_forces = [f for f in omm_system.getForces() if type(f) == mm.CustomNonbondedForce]
        pair12, pair13, pair14 = system._topology.get_12_13_14_pairs()
        for atom1, atom2 in pair12 + pair13:
            nbforce.addException(atom1.id, atom2.id, 0.0, 1.0, 0.0)
            for f in custom_nb_forces:
                f.addExclusion(atom1.id, atom2.id)
        # As long as 1-4 LJ OR Coulomb need to be scaled, then this pair should be excluded from ALL non-bonded forces.
        # This is required by OpenMM's internal implementation.
        # Even though NonbondedForce can handle 1-4 vdW, we use it only for 1-4 Coulomb.
        # And use CustomBondForce to handle 1-4 vdW, which makes it more clear for energy decomposition.
        if system._ff.scale_14_vdw != 1 or system._ff.scale_14_coulomb != 1:
            pair14_forces = {}  # {VdwTerm: mm.NbForce}
            for atom1, atom2 in pair14:
                charge_prod = atom1.charge * atom2.charge * system._ff.scale_14_coulomb
                nbforce.addException(atom1.id, atom2.id, charge_prod, 1.0, 0.0)
                for f in custom_nb_forces:
                    f.addExclusion(atom1.id, atom2.id)
                if system._ff.scale_14_vdw == 0:
                    continue
                vdw = system._ff.get_vdw_term(atom1.type, atom2.type)
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
                        cbforce.setUsesPeriodicBoundaryConditions(True)
                        cbforce.setForceGroup(ForceGroup.VDW)
                        omm_system.addForce(cbforce)
                        pair14_forces[MieTerm] = cbforce
                    epsilon = vdw.epsilon * system._ff.scale_14_vdw
                    if type(vdw) == LJ126Term:
                        cbforce.addBond(atom1.id, atom2.id, [epsilon, vdw.sigma, 12, 6])
                    elif type(vdw) == MieTerm:
                        cbforce.addBond(atom1.id, atom2.id,
                                        [epsilon, vdw.sigma, vdw.repulsion, vdw.attraction])
                else:
                    raise Exception('1-4 scaling for vdW terms other than LJ126Term and MieTerm '
                                    'haven\'t been implemented')

        ### Set up Drude particles ##############################################################
        for polar_class in system.polarizable_classes:
            if polar_class == DrudeTerm:
                logger.info('Setting up Drude polarizations...')
                pforce = mm.DrudeForce()
                pforce.setForceGroup(ForceGroup.DRUDE)
                omm_system.addForce(pforce)
                parent_idx_thole = {}  # {parent: (index in DrudeForce, thole)} for addScreenPair
                for parent, drude in system.drude_pairs.items():
                    pterm = system.polarizable_terms[parent]
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
                drude_exclusions14 = []
                for atom1, atom2 in pair14:
                    drude1 = system.drude_pairs.get(atom1)
                    drude2 = system.drude_pairs.get(atom2)
                    if drude1 is not None:
                        drude_exclusions14.append((drude1, atom2))
                    if drude2 is not None:
                        drude_exclusions14.append((atom1, drude2))
                    if drude1 is not None and drude2 is not None:
                        drude_exclusions14.append((drude1, drude2))
                for a1, a2 in drude_exclusions14:
                    charge_prod = a1.charge * a2.charge * system._ff.scale_14_coulomb
                    nbforce.addException(a1.id, a2.id, charge_prod, 1.0, 0.0)
                    for f in custom_nb_forces:
                        f.addExclusion(a1.id, a2.id)
            else:
                raise Exception('Polarizable terms other that DrudeTerm haven\'t been implemented')

        ### Remove COM motion ###################################################################
        logger.info('Setting up COM motion remover...')
        omm_system.addForce(mm.CMMotionRemover(10))

        return omm_system
