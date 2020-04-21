import itertools
import numpy as np
import warnings
from ..topology import Topology, Atom, UnitCell, Psf, Bond, Angle, Dihedral, Improper
from ..trajectory import Frame, Trajectory, Gro
from ..forcefield import FFSet
from ..forcefield.ffterm import *
from ..constant import *

try:
    from simtk import openmm as mm
    from simtk.openmm import app
except ImportError:
    OMM_EXIST = False
else:
    OMM_EXIST = True


class System():
    def __init__(self, topology: Topology, params: FFSet, positions=None, cell=None, verbose=False):
        self._topology = topology
        self._params = params
        self._verbose = verbose
        self._frame = Frame(topology.n_atom)

        if positions is not None:
            self._frame.positions = positions
        elif topology.has_position:
            self._frame.positions = topology.positions
        else:
            raise Exception('Positions should be provided with topology or positions')

        if cell is not None:
            self._frame.cell = cell
        else:
            self._frame.cell = topology.cell

        self.extract_required_ff_terms()

    def extract_required_ff_terms(self):
        self.atom_types: [AtomType] = []
        self.type_names: [str] = []
        self.vdw_terms: [VdwTerm] = []  # all required self-self and pairwise vdW terms
        self.bond_terms: {int: BondTerm} = {}  # key is id(Bond)
        self.angle_terms: {int: AngleTerm} = {}  # key is id(Angle)
        self.dihedral_terms: {int: DihedralTerm} = {}  # key is id(Dihedral)
        self.improper_terms: {int: ImproperTerm} = {}  # key is id(Improper)
        self.polarizable_terms: {str: PolarizableTerm} = {}  # key is Atom.type
        self.drude_pairs: {Atom: Atom} = dict(self._topology.get_drude_pairs())  # {parent: drude}
        self.constrain_bonds: {int: float} = {}  # key is id(Bond), value is distance
        self.constrain_angles: {int: float} = {}  # key is id(Angle), value is 1-3 distance

        for atom in self._topology.atoms:
            atype = self._params.atom_types[atom.type]
            if atype is None:
                raise Exception(f'AtomType {atom.type} not found in FF')

            if atype not in self.atom_types:
                self.atom_types.append(atype)
                self.type_names.append(atype.name)
                self.vdw_terms.append(self._params.get_vdw_term(atype, atype))

        for atype1, atype2 in itertools.combinations(self.atom_types, 2):
            vdw = self._params.pairwise_vdw_terms.get(VdwTerm(atype1.eqt_vdw, atype2.eqt_vdw).name)
            if vdw is not None:
                self.vdw_terms.append(vdw)

        for bond in self._topology.bonds:
            if bond.is_drude:
                # the Drude-parent bonds will be handled by DrudeForce below
                continue
            at1 = self._params.atom_types[bond.atom1.type].eqt_bond
            at2 = self._params.atom_types[bond.atom2.type].eqt_bond
            bterm = BondTerm(at1, at2, 0)
            try:
                bterm = self._params.bond_terms[bterm.name]
            except:
                raise Exception(f'{str(bterm)} for {str(bond)} not found in FF')
            self.bond_terms[id(bond)] = bterm
            if bterm.fixed:
                self.constrain_bonds[id(bond)] = bterm.length

        for angle in self._topology.angles:
            at1 = self._params.atom_types[angle.atom1.type].eqt_ang_s
            at2 = self._params.atom_types[angle.atom2.type].eqt_ang_c
            at3 = self._params.atom_types[angle.atom3.type].eqt_ang_s
            aterm = AngleTerm(at1, at2, at3, 0)
            try:
                aterm = self._params.angle_terms[aterm.name]
            except:
                raise Exception(f'{str(aterm)} for {str(angle)} not found in FF')
            self.angle_terms[id(angle)] = aterm
            if aterm.fixed:
                bond1 = next(b for b in self._topology.bonds if b == Bond(angle.atom1, angle.atom2))
                bond2 = next(b for b in self._topology.bonds if b == Bond(angle.atom2, angle.atom3))
                # constrain the angle only if two bonds are also constrained
                if id(bond1) in self.constrain_bonds and id(bond2) in self.constrain_bonds:
                    d1, d2 = self.constrain_bonds[id(bond1)], self.constrain_bonds[id(bond2)]
                    self.constrain_angles[id(angle)] = math.sqrt(
                        d1 * d1 + d2 * d2 - 2 * d1 * d2 * math.cos(aterm.theta * PI / 180))

        for dihedral in self._topology.dihedrals:
            at1 = self._params.atom_types[dihedral.atom1.type].eqt_dih_s
            at2 = self._params.atom_types[dihedral.atom2.type].eqt_dih_c
            at3 = self._params.atom_types[dihedral.atom3.type].eqt_dih_c
            at4 = self._params.atom_types[dihedral.atom4.type].eqt_dih_s
            for a1, a2, a3, a4 in [(at1, at2, at3, at4),
                                   (at1, at2, at3, '*'),
                                   ('*', at2, at3, at4),
                                   ('*', at2, at3, '*')]:
                dterm = self._params.dihedral_terms.get(DihedralTerm(a1, a2, a3, a4).name)
                if dterm is not None:
                    self.dihedral_terms[id(dihedral)] = dterm
                    break
            else:
                dterm = DihedralTerm(at1, at2, at3, at4)
                raise Exception(f'{str(dterm)} for {str(dihedral)} not found in FF')

        for improper in self._topology.impropers:
            at1 = self._params.atom_types[improper.atom1.type].eqt_imp_c
            at2 = self._params.atom_types[improper.atom2.type].eqt_imp_s
            at3 = self._params.atom_types[improper.atom3.type].eqt_imp_s
            at4 = self._params.atom_types[improper.atom4.type].eqt_imp_s
            for a2, a3, a4 in [(at2, at3, at4),
                               (at2, at3, '*'),
                               (at2, at4, '*'),
                               (at3, at4, '*'),
                               (at2, '*', '*'),
                               (at3, '*', '*'),
                               (at4, '*', '*'),
                               ('*', '*', '*')]:
                iterm = self._params.improper_terms.get(ImproperTerm(at1, a2, a3, a4).name)
                if iterm is not None:
                    self.improper_terms[id(improper)] = iterm
                    break
            else:
                iterm = ImproperTerm(at1, at2, at3, at4)
                raise Exception(f'{str(iterm)} for {str(improper)} not found in FF')

        for parent in self.drude_pairs.keys():
            pterm = PolarizableTerm(self._params.atom_types[parent.type].eqt_polar)
            try:
                self.polarizable_terms[parent.type] = self._params.polarizable_terms[pterm.name]
            except:
                raise Exception(f'{str(pterm)} for {str(parent)} not found')

        self.vdw_classes = set([term.__class__ for term in self.vdw_terms])
        self.bond_classes = set([term.__class__ for term in self.bond_terms.values()])
        self.angle_classes = set([term.__class__ for term in self.angle_terms.values()])
        self.dihedral_classes = set([term.__class__ for term in self.dihedral_terms.values()])
        self.improper_classes = set([term.__class__ for term in self.improper_terms.values()])
        self.polarizable_classes = set([term.__class__ for term in self.polarizable_terms.values()])

        if self._verbose:
            print('%10i atom types' % len(self.atom_types))
            print('%10i bonds constrained' % len(self.constrain_bonds))
            print('%10i angles constrained' % len(self.constrain_angles))

    def export_lmp(self):
        pass

    def export_gmx(self):
        pass

    def to_omm_system(self):
        if not OMM_EXIST:
            raise ImportError('Can not import OpenMM')

        omm_system = mm.System()
        omm_system.setDefaultPeriodicBoxVectors(*self._frame.cell.vectors)
        for atom in self._topology.atoms:
            omm_system.addParticle(atom.mass)

        ### Set up bonds #######################################################################
        if self.bond_classes - {HarmonicBondTerm} > set():
            raise Exception('Bond terms other that HarmonicBondTerm haven\'t been implemented')
        bforce = mm.HarmonicBondForce()
        bforce.setUsesPeriodicBoundaryConditions(True)
        bforce.setForceGroup(1)
        omm_system.addForce(bforce)
        for bond in self._topology.bonds:
            if bond.is_drude:
                continue
            bterm = self.bond_terms[id(bond)]
            bforce.addBond(bond.atom1.id, bond.atom2.id, bterm.length, bterm.k * 2)

        ### Set up angles #######################################################################
        for angle_class in self.angle_classes:
            if angle_class == HarmonicAngleTerm:
                aforce = mm.HarmonicAngleForce()
                for angle in self._topology.angles:
                    aterm = self.angle_terms[id(angle)]
                    if type(aterm) == HarmonicAngleTerm:
                        aforce.addAngle(angle.atom1.id, angle.atom2.id, angle.atom3.id,
                                        aterm.theta * PI / 180, aterm.k * 2)
            elif angle_class == SDKAngleTerm:
                aforce = mm.CustomCompoundBondForce(
                    'k*(theta-theta0)^2+step(rmin-r)*LJ96;'
                    'LJ96=6.75*epsilon*((sigma/r)^9-(sigma/r)^6)+epsilon;'
                    'theta=angle(p1,p2,p3);'
                    'r=distance(p1,p3);'
                    'rmin=1.144714*sigma')
                aforce.addPerBondParameter('theta0')
                aforce.addPerBondParameter('k')
                aforce.addPerBondParameter('epsilon')
                aforce.addPerBondParameter('sigma')
                for angle in self._topology.angles:
                    aterm = self.angle_terms[id(angle)]
                    if type(aterm) == SDKAngleTerm:
                        aforce.addAngle(
                            [angle.atom1.id, angle.atom2.id, angle.atom3.id],
                            [aterm.theta * PI / 180, aterm.k, aterm.epsilon, aterm.sigma])
            else:
                raise Exception('Angle terms other that HarmonicAngleTerm and SDKAngleTerm '
                                'haven\'t been implemented')
            aforce.setUsesPeriodicBoundaryConditions(True)
            aforce.setForceGroup(2)
            omm_system.addForce(aforce)

        ### Set up constraints #################################################################
        for bond in self._topology.bonds:
            if id(bond) in self.constrain_bonds:
                omm_system.addConstraint(bond.atom1.id, bond.atom2.id,
                                         self.constrain_bonds[id(bond)])
        for angle in self._topology.angles:
            if id(angle) in self.constrain_angles:
                omm_system.addConstraint(angle.atom1.id, angle.atom3.id,
                                         self.constrain_angles[id(angle)])

        ### Set up dihedrals ###################################################################
        dforce = mm.PeriodicTorsionForce()
        dforce.setUsesPeriodicBoundaryConditions(True)
        dforce.setForceGroup(3)
        omm_system.addForce(dforce)
        for dihedral in self._topology.dihedrals:
            dterm = self.dihedral_terms[id(dihedral)]
            ia1, ia2, ia3, ia4 = dihedral.atom1.id, dihedral.atom2.id, dihedral.atom3.id, dihedral.atom4.id
            if type(dterm) == PeriodicDihedralTerm:
                for par in dterm.parameters:
                    dforce.addTorsion(ia1, ia2, ia3, ia4, par.n, par.phi * PI / 180, par.k)
            elif type(dterm) == OplsDihedralTerm:
                if dterm.k1 != 0:
                    dforce.addTorsion(ia1, ia2, ia3, ia4, 1, 0, dterm.k1)
                if dterm.k2 != 0:
                    dforce.addTorsion(ia1, ia2, ia3, ia4, 2, PI, dterm.k2)
                if dterm.k3 != 0:
                    dforce.addTorsion(ia1, ia2, ia3, ia4, 3, 0, dterm.k3)
                if dterm.k4 != 0:
                    dforce.addTorsion(ia1, ia2, ia3, ia4, 4, PI, dterm.k4)
            else:
                raise Exception(
                    'Dihedral terms other that PeriodicDihedralTerm and OplsDihedralTerm'
                    'haven\'t been implemented')

        ### Set up impropers ####################################################################
        for improper_class in self.improper_classes:
            if improper_class == PeriodicImproperTerm:
                iforce = mm.CustomTorsionForce('k*(1+cos(2*theta-phi0))')
                iforce.addPerTorsionParameter('phi0')
                iforce.addPerTorsionParameter('k')
                for improper in self._topology.impropers:
                    iterm = self.improper_terms[id(improper)]
                    if type(iterm) == PeriodicImproperTerm:
                        iforce.addTorsion(improper.atom1.id, improper.atom2.id,
                                          improper.atom3.id, improper.atom4.id,
                                          [iterm.phi * PI / 180, iterm.k])
            elif improper_class == HarmonicImproperTerm:
                iforce = mm.CustomTorsionForce(f'k*min(dtheta,2*pi-dtheta)^2;'
                                               f'dtheta=abs(theta-phi0);'
                                               f'pi={PI}')
                iforce.addPerTorsionParameter('phi0')
                iforce.addPerTorsionParameter('k')
                for improper in self._topology.impropers:
                    iterm = self.improper_terms[id(improper)]
                    if type(iterm) == HarmonicImproperTerm:
                        iforce.addTorsion(improper.atom1.id, improper.atom2.id,
                                          improper.atom3.id, improper.atom4.id,
                                          [iterm.phi * PI / 180, iterm.k])
            else:
                raise Exception(
                    'Improper terms other that PeriodicImproperTerm and HarmonicImproperTerm '
                    'haven\'t been implemented')
            iforce.setUsesPeriodicBoundaryConditions(True)
            iforce.setForceGroup(4)
            omm_system.addForce(iforce)

        ### Set up coulomb interactions #########################################################
        n_type = len(self.atom_types)
        cutoff = self._params.vdw_cutoff
        # NonbonedForce is not flexible enough. Use it only for Coulomb interactions
        # CustomNonbondedForce handles vdW interactions
        nbforce = mm.NonbondedForce()
        nbforce.setNonbondedMethod(mm.NonbondedForce.PME)
        nbforce.setCutoffDistance(cutoff)
        nbforce.setUseDispersionCorrection(False)
        try:
            nbforce.setExceptionsUsePeriodicBoundaryConditions(True)
        except:
            warnings.warn('Cannot apply PBC for vdW exceptions. '
                          'Be careful if there\'s infinite structures in the system')
        nbforce.setForceGroup(5)
        omm_system.addForce(nbforce)
        for atom in self._topology.atoms:
            nbforce.addParticle(atom.charge, 1.0, 0.0)

        ### Set up vdW interactions #########################################################
        for vdw_class in self.vdw_classes:
            if vdw_class == LJ126Term:
                if self._params.vdw_long_range == FFSet.VDW_LONGRANGE_SHIFT:
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
                for idx in range(n_type):
                    for j in range(n_type):
                        at1 = self.atom_types[idx]
                        at2 = self.atom_types[j]
                        vdw = self._params.get_vdw_term(at1, at2)
                        if type(vdw) == LJ126Term:
                            A = 4 * vdw.epsilon * vdw.sigma ** 12
                            B = 4 * vdw.epsilon * vdw.sigma ** 6
                        else:
                            A = B = 0
                        A_list[idx + n_type * j] = A
                        B_list[idx + n_type * j] = B
                cforce.addTabulatedFunction('A', mm.Discrete2DFunction(n_type, n_type, A_list))
                cforce.addTabulatedFunction('B', mm.Discrete2DFunction(n_type, n_type, B_list))

                for atom in self._topology.atoms:
                    id_type = self.type_names.index(atom.type)
                    cforce.addParticle([id_type])

            elif vdw_class == MieTerm:
                if self._params.vdw_long_range == FFSet.VDW_LONGRANGE_SHIFT:
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
                for idx in range(n_type):
                    for j in range(n_type):
                        at1 = self.atom_types[idx]
                        at2 = self.atom_types[j]
                        vdw = self._params.get_vdw_term(at1, at2)
                        if type(vdw) == MieTerm:
                            A = vdw.factor_energy() * vdw.epsilon * vdw.sigma ** vdw.repulsion
                            B = vdw.factor_energy() * vdw.epsilon * vdw.sigma ** vdw.attraction
                            REP = vdw.repulsion
                            ATT = vdw.attraction
                            SHIFT = A / cutoff ** REP - B / cutoff ** ATT
                        else:
                            A = B = REP = ATT = SHIFT = 0
                        A_list[idx + n_type * j] = A
                        B_list[idx + n_type * j] = B
                        REP_list[idx + n_type * j] = REP
                        ATT_list[idx + n_type * j] = ATT
                        SHIFT_list[idx + n_type * j] = SHIFT
                cforce.addTabulatedFunction('A', mm.Discrete2DFunction(n_type, n_type, A_list))
                cforce.addTabulatedFunction('B', mm.Discrete2DFunction(n_type, n_type, B_list))
                cforce.addTabulatedFunction('REP',
                                            mm.Discrete2DFunction(n_type, n_type, REP_list))
                cforce.addTabulatedFunction('ATT',
                                            mm.Discrete2DFunction(n_type, n_type, ATT_list))
                if self._params.vdw_long_range == FFSet.VDW_LONGRANGE_SHIFT:
                    cforce.addTabulatedFunction('SHIFT',
                                                mm.Discrete2DFunction(n_type, n_type, SHIFT_list))

                for atom in self._topology.atoms:
                    id_type = self.type_names.index(atom.type)
                    cforce.addParticle([id_type])

            else:
                raise Exception('vdW terms other than LJ126Term and MieTerm '
                                'haven\'t been implemented')
            cforce.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            cforce.setCutoffDistance(cutoff)
            if self._params.vdw_long_range == FFSet.VDW_LONGRANGE_CORRECT:
                cforce.setUseLongRangeCorrection(True)
            else:
                cforce.setUseLongRangeCorrection(False)
            cforce.setForceGroup(6)
            omm_system.addForce(cforce)

        ### Set up 1-2, 1-3 and 1-4 exceptions ##################################################
        custom_nb_forces = [f for f in omm_system.getForces() if type(f) == mm.CustomNonbondedForce]
        pair12, pair13, pair14 = self._topology.get_12_13_14_pairs()
        for atom1, atom2 in pair12 + pair13:
            nbforce.addException(atom1.id, atom2.id, 0.0, 1.0, 0.0)
            for f in custom_nb_forces:
                f.addExclusion(atom1.id, atom2.id)
        # As long as 1-4 LJ OR Coulomb need to be scaled, then this pair should be excluded from ALL non-bonded forces.
        # This is required by OpenMM's internal implementation.
        # Even though NonbondedForce can handle 1-4 vdW, we use it only for 1-4 Coulomb.
        # And use CustomBondForce to handle 1-4 vdW, which makes it more clear for energy decomposition.
        if self._params.scale_14_vdw != 1 or self._params.scale_14_coulomb != 1:
            pair14_forces = {}  # {VdwTerm: mm.NbForce}
            for atom1, atom2 in pair14:
                charge_prod = atom1.charge * atom2.charge * self._params.scale_14_coulomb
                nbforce.addException(atom1.id, atom2.id, charge_prod, 1.0, 0.0)
                for f in custom_nb_forces:
                    f.addExclusion(atom1.id, atom2.id)
                if self._params.scale_14_vdw == 0:
                    continue
                vdw = self._params.get_vdw_term(self._params.atom_types[atom1.type],
                                                self._params.atom_types[atom2.type])
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
                        cbforce.setForceGroup(6)
                        omm_system.addForce(cbforce)
                        pair14_forces[MieTerm] = cbforce
                    epsilon = vdw.epsilon * self._params.scale_14_vdw
                    if type(vdw) == LJ126Term:
                        cbforce.addBond(atom1.id, atom2.id, [epsilon, vdw.sigma, 12, 6])
                    elif type(vdw) == MieTerm:
                        cbforce.addBond(atom1.id, atom2.id,
                                        [epsilon, vdw.sigma, vdw.repulsion, vdw.attraction])
                else:
                    raise Exception('1-4 scaling for vdW terms other than LJ126Term and MieTerm '
                                    'haven\'t been implemented')

        ### Set up Drude particles ##############################################################
        if self.polarizable_classes - {DrudeTerm} > set():
            raise Exception('Polarizable terms other that DrudePolarizableTerm '
                            'haven\'t been implemented')
        pforce = mm.DrudeForce()
        pforce.setForceGroup(7)
        omm_system.addForce(pforce)
        parent_idx_thole = {}  # {parent: (index in DrudeForce, thole)} for addScreenPair
        for parent, drude in self.drude_pairs.items():
            pterm = self.polarizable_terms[parent.type]
            n_H = len([atom for atom in parent.bond_partners if atom.symbol == 'H'])
            alpha = pterm.alpha + n_H * pterm.merge_alpha_H
            idx = pforce.addParticle(drude.id, parent.id, -1, -1, -1, drude.charge, alpha, 0, 0)
            parent_idx_thole[parent] = (idx, pterm.thole)

        # exclude the non-boned interactions between Drude and parent
        # and those concerning Drude particles in 1-2 and 1-3 pairs
        # pairs formed by real atoms have already been handled above
        # also apply thole screening between 1-2 and 1-3 Drude dipole pairs
        drude_exclusions = list(self.drude_pairs.items())
        for atom1, atom2 in pair12 + pair13:
            drude1 = self.drude_pairs.get(atom1)
            drude2 = self.drude_pairs.get(atom2)
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
            drude1 = self.drude_pairs.get(atom1)
            drude2 = self.drude_pairs.get(atom2)
            if drude1 is not None:
                drude_exclusions14.append((drude1, atom2))
            if drude2 is not None:
                drude_exclusions14.append((atom1, drude2))
            if drude1 is not None and drude2 is not None:
                drude_exclusions14.append((drude1, drude2))
        for a1, a2 in drude_exclusions14:
            charge_prod = a1.charge * a2.charge * self._params.scale_14_coulomb
            nbforce.addException(a1.id, a2.id, charge_prod, 1.0, 0.0)
            for f in custom_nb_forces:
                f.addExclusion(a1.id, a2.id)

        return omm_system
