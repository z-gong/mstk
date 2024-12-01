import itertools
from .system import System
from ..forcefield import *
from ..topology import *
from ..trajectory import *
from .. import logger


class GromacsExporter:
    '''
    GromacsExporter export a :class:`System` to input files for Gromacs.

    The following potential functions are currently supported:

    * :class:`~mstk.forcefield.LJ126Term`
    * :class:`~mstk.forcefield.HarmonicBondTerm`
    * :class:`~mstk.forcefield.MorseBondTerm`
    * :class:`~mstk.forcefield.HarmonicAngleTerm`
    * :class:`~mstk.forcefield.HarmonicCosineAngleTerm`
    * :class:`~mstk.forcefield.LinearAngleTerm`
    * :class:`~mstk.forcefield.OplsDihedralTerm`
    * :class:`~mstk.forcefield.PeriodicDihedralTerm`
    * :class:`~mstk.forcefield.OplsImproperTerm`
    * :class:`~mstk.forcefield.HarmonicImproperTerm`
    * :class:`~mstk.forcefield.DrudePolarTerm`

    Note that the Drude polarization in GROMACS have never been thoroughly tested. Use it with care.

    In order to use GROMACS for trajectory analysis, the following potential functions can also be exported.
    But the exported files are only used to generate tpr file, which is required by most gmx analysis commands.
    Do NOT use the exported files for production simulations.

    * :class:`~mstk.forcefield.MieTerm` will be exported in the LJ-12-6 form
    * :class:`~mstk.forcefield.MorseVdwTerm` will be exported in the LJ-12-6 form
    * :class:`~mstk.forcefield.QuarticBondTerm` will be exported in the harmonic form
    * :class:`~mstk.forcefield.QuarticAngleTerm` will be exported in the harmonic form
    * :class:`~mstk.forcefield.SDKAngleTerm` will be exported in the harmonic form
    '''

    def __init__(self):
        pass

    @staticmethod
    def export(system, gro_out, top_out, mdp_out, **kwargs):
        '''
        Generate input files for Gromacs from a system

        # TODO Export virtual sites other than TIP4P

        Parameters
        ----------
        system : System
        gro_out : str or None
        top_out : str or None
        mdp_out : str or None
        '''
        if not system.use_pbc:
            raise Exception('PBC required for exporting GROMACS')

        supported_terms = {LJ126Term, MieTerm, MorseVdwTerm,
                           HarmonicBondTerm, QuarticBondTerm, MorseBondTerm,
                           HarmonicAngleTerm, QuarticAngleTerm, HarmonicCosineAngleTerm, SDKAngleTerm, LinearAngleTerm,
                           OplsDihedralTerm, PeriodicDihedralTerm,
                           OplsImproperTerm, HarmonicImproperTerm,
                           DrudePolarTerm}
        unsupported = system.ff.energy_term_classes - supported_terms
        if unsupported != set():
            raise Exception('Unsupported FF terms: %s' % (', '.join(map(lambda x: x.__name__, unsupported))))

        if MieTerm in system.ff.vdw_term_classes:
            logger.warning('MieTerm not supported by GROMACS. Exported in LJ-12-6 form')
        if MorseVdwTerm in system.ff.vdw_term_classes:
            logger.warning('MorseVdwTerm not supported by GROMACS. Exported in LJ-12-6 form')
        if QuarticBondTerm in system.ff.bond_term_classes:
            logger.warning('QuarticBondTerm not supported by GROMACS. Exported in harmonic form')
        if QuarticAngleTerm in system.ff.angle_term_classes:
            logger.warning('QuarticAngleTerm not supported by GROMACS. Exported in harmonic form')
        if SDKAngleTerm in system.ff.angle_term_classes:
            logger.warning('SDKAngleTerm not supported by GROMACS. Exported in harmonic form')
        if DrudePolarTerm in system.ff.polar_term_classes:
            logger.warning('DrudePolarTerm not well tested in GROMACS. Use it with caution')

        if system.topology.virtual_site_classes - {TIP4PSite} != set():
            raise Exception('Virtual sites other than TIP4PSite haven\'t been implemented')

        if gro_out is not None:
            GromacsExporter._export_gro(system, gro_out)
        if top_out is not None:
            GromacsExporter._export_top(system, top_out)
        if mdp_out is not None:
            GromacsExporter._export_mdp(system, mdp_out)

    @staticmethod
    def _export_gro(system: System, gro_out='conf.gro'):
        frame = Frame(system._topology.n_atom)
        frame.cell = system._topology.cell
        frame.positions = system._topology.positions
        gro = Trajectory(gro_out, 'w')
        gro.write_frame(frame, system._topology)
        gro.close()

    @staticmethod
    def _export_top(system: System, top_out='topol.top'):
        mols_unique = system._topology.get_unique_molecules(deepcopy=False)

        string = '; GROMACS topol file created by mstk\n'

        string += '\n[ defaults ]\n'
        string += '; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n'
        _comb = 2 if system.ff.lj_mixing_rule == ForceField.LJ_MIXING_LB else 3
        string += '%6i %6i %12s %10.4f %10.4f\n' % (
            1, _comb, 'yes', system.ff.scale_14_vdw, system.ff.scale_14_coulomb)

        string += '\n[ atomtypes ]\n'
        string += ';     name       mass     charge      ptype      sigma      epsilon\n'
        drude_types = set()
        if DrudePolarTerm in system.ff.polar_term_classes:
            for atom in system.topology.atoms:
                if atom.is_drude:
                    drude_types.add(atom.type)
        for atype in system.ff.atom_types.values():
            ptype = 'S' if atype.name in drude_types else 'A'
            vdw = system.ff.get_vdw_term(atype, atype)
            if vdw.__class__ in (LJ126Term, MieTerm):
                string += '%10s %10.4f %12.6f %6s %12.6f %12.6f  ; %s\n' % (
                    atype.name, 0.0, 0.0, ptype, vdw.sigma, vdw.epsilon, ' '.join(vdw.comments))
            elif vdw.__class__ == MorseVdwTerm:
                string += '%10s %10.4f %12.6f %6s %12.6f %12.6f  ; %s\n' % (
                    atype.name, 0.0, 0.0, ptype, vdw.r0, vdw.depth, ' '.join(vdw.comments))
            else:
                raise Exception('Unsupported vdW term')

        string += '\n[ nonbond_params ]\n'
        string += ';       i        j        func      sigma      epsilon\n'
        for at1, at2 in itertools.combinations(system.ff.atom_types.values(), 2):
            try:
                vdw = system.ff.get_vdw_term(at1, at2, mixing=False)
            except:
                continue
            if vdw.__class__ in (LJ126Term, MieTerm):
                string += '%10s %10s %6i %12.6f %12.6f  ; %s\n' % (
                    at1.name, at2.name, 1, vdw.sigma, vdw.epsilon, ' '.join(vdw.comments))
            elif vdw.__class__ == MorseVdwTerm:
                string += '%10s %10s %6i %12.6f %12.6f  ; %s\n' % (
                    at1.name, at2.name, 1, vdw.r0, vdw.depth, ' '.join(vdw.comments))
            else:
                raise Exception('Unsupported vdW term')

        for i, mol in enumerate(mols_unique.keys()):
            mol: Molecule
            string += '\n[ moleculetype ]\n'
            string += '; Name            nrexcl\n'
            string += f'{mol.name}_{i + 1}     3\n'

            string += '\n[ atoms ]\n'
            string += ';%5s %10s %6s %10s %6s %6s %12s %10s\n' % (
                'nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass')
            for j, atom in enumerate(mol.atoms):
                # GROMACS use SCF approach for Drude simulation. Set mass of Drude to zero
                mass = 0 if atom.is_drude else atom.mass
                if atom in system.drude_pairs:
                    mass += system.drude_pairs[atom].mass
                string += '%6i %10s %6i %10s %6s %6i %12.6f %10.4f\n' % (
                    j + 1, atom.type, atom.residue.id_in_mol + 1, atom.residue.name,
                    atom.symbol, j + 1, atom.charge, mass)

            if TIP4PSite in system.topology.virtual_site_classes:
                string += '\n[ virtual_sites3 ]\n'
                string += ';%5s %6s %6s %6s %6s %12s %12s\n' % ('Vsite', 'from1', 'from2', 'from3', 'funct', 'a', 'b')
                for atom in mol.atoms:
                    vsite = atom.virtual_site
                    if type(vsite) is TIP4PSite:
                        O, H1, H2 = vsite.parents
                        w_O, w_H1, w_H2 = system.get_TIP4P_linear_coeffs(atom)
                        string += '%6i %6i %6i %6i %6i %12.6f %12.6f\n' % (
                            atom.id_in_mol + 1, O.id_in_mol + 1, H1.id_in_mol + 1, H2.id_in_mol + 1, 1, w_H1, w_H2)
                    elif vsite is not None:
                        raise Exception('Virtual sites other than TIP4PSite haven\'t been implemented')

            string += '\n[ pairs ]\n'
            pairs12, pairs13, pairs14 = mol.get_12_13_14_pairs()
            for a1, a2 in pairs14:
                string += '%6i %6i %6i\n' % (a1.id_in_mol + 1, a2.id_in_mol + 1, 1)

            # 1-4 pairs concerning Drude particles
            drude_pairs = dict(mol.get_drude_pairs())
            _pairs14_drude = []
            for a1, a2 in pairs14:
                if a1 in drude_pairs:
                    _pairs14_drude.append((drude_pairs[a1], a2))
                if a2 in drude_pairs:
                    _pairs14_drude.append((a1, drude_pairs[a2]))
                if a1 in drude_pairs and a2 in drude_pairs:
                    _pairs14_drude.append((drude_pairs[a1], drude_pairs[a2]))

            # 1-4 pairs concerning virtual sites
            vsite_pairs = dict(mol.get_virtual_site_pairs())
            _pairs14_vsite = []
            for atom1, atom2 in pairs14:
                vsite1 = vsite_pairs.get(atom1)
                vsite2 = vsite_pairs.get(atom2)
                drude1 = drude_pairs.get(atom1)
                drude2 = drude_pairs.get(atom2)
                if vsite1 is not None:
                    _pairs14_vsite.append((vsite1, atom2))
                    if drude2 is not None:
                        _pairs14_vsite.append((vsite1, drude2))
                if vsite2 is not None:
                    _pairs14_vsite.append((vsite2, atom1))
                    if drude1 is not None:
                        _pairs14_vsite.append((vsite2, drude1))
                if None not in [vsite1, vsite2]:
                    _pairs14_vsite.append((vsite1, vsite2))

            for a1, a2 in _pairs14_drude + _pairs14_vsite:
                string += '%6i %6i %6i\n' % (a1.id_in_mol + 1, a2.id_in_mol + 1, 1)

            ### Drude polarization #############################################
            if DrudePolarTerm in system.ff.polar_term_classes:
                string += '\n[ polarization ]\n'
                string += ';  ai    aj   func    alpha     delta     khyp\n'
                for parent, drude in drude_pairs.items():
                    string += '%5d %5d   %4d  %9.6f  %7.4f  %9.4e\n' % (
                        parent.id_in_mol + 1, drude.id_in_mol + 1,
                        2, parent.alpha, 0.02, 16.768e8)

                string += '\n[ thole_polarization ]\n'
                string += ';  ai    aj    ak    al   func   thole   alpha1   alpha2\n'
                for a1, a2 in pairs12 + pairs13:
                    if a1 in drude_pairs and a2 in drude_pairs:
                        string += '%5i %5i %5i %5i   %4i   %9.5f  %9.6f  %9.6f\n' % (
                            a1.id_in_mol + 1, drude_pairs[a1].id_in_mol + 1,
                            a2.id_in_mol + 1, drude_pairs[a2].id_in_mol + 1,
                            1, (a1.thole + a2.thole) / 2, a1.alpha, a2.alpha)

            string += '\n[ constraints ]\n'
            for bond in mol.bonds:
                distance = system.constrain_bonds.get(bond)
                if distance is not None:
                    a1, a2 = bond.atom1, bond.atom2
                    string += '%6i %6i %6i %12.6f\n' % (
                        a1.id_in_mol + 1, a2.id_in_mol + 1, 1, distance)
            for angle in mol.angles:
                distance = system.constrain_angles.get(angle)
                if distance is not None:
                    a1, a3 = angle.atom1, angle.atom3
                    string += '%6i %6i %6i %12.6f\n' % (
                        a1.id_in_mol + 1, a3.id_in_mol + 1, 2, distance)

            string += '\n[ bonds ]\n'
            for bond in mol.bonds:
                if bond.is_drude:
                    continue
                if bond in system.constrain_bonds:
                    continue
                bterm = system.bond_terms[bond]
                if bterm.__class__ == HarmonicBondTerm:
                    a1, a2 = bond.atom1, bond.atom2
                    string += '%6i %6i %6i %12.6f %12.4f\n' % (
                        a1.id_in_mol + 1, a2.id_in_mol + 1,
                        1, bterm.length, bterm.k * 2)
                elif bterm.__class__ == QuarticBondTerm:
                    a1, a2 = bond.atom1, bond.atom2
                    string += '%6i %6i %6i %12.6f %12.4f\n' % (
                        a1.id_in_mol + 1, a2.id_in_mol + 1,
                        1, bterm.length, bterm.k2 * 2)
                elif bterm.__class__ == MorseBondTerm:
                    a1, a2 = bond.atom1, bond.atom2
                    string += '%6i %6i %6i %12.6f %12.4f %12.4f\n' % (
                        a1.id_in_mol + 1, a2.id_in_mol + 1,
                        3, bterm.length, bterm.depth, (bterm.k / bterm.depth) ** 0.5)
                else:
                    raise Exception('Unsupported bond term')

            # Exclusions concerning Drude particles
            _exclusions_drude = {}  # {drude1: [atom2, drude2, atom3, drude3, ...]}
            for a1, a2 in pairs12 + pairs13 + pairs14:
                if a1 in drude_pairs:
                    d1 = drude_pairs[a1]
                    if d1 not in _exclusions_drude:
                        _exclusions_drude[d1] = []
                    _exclusions_drude[d1].append(a2)
                if a2 in drude_pairs:
                    d2 = drude_pairs[a2]
                    if d2 not in _exclusions_drude:
                        _exclusions_drude[d2] = []
                    _exclusions_drude[d2].append(a1)
                if a1 in drude_pairs and a2 in drude_pairs:
                    _exclusions_drude[d1].append(d2)

            # Exclusions concerning virtual sites
            _exclusions_vsite = {}  # {avsite1: [atom1, drude1, atom2, drude2, avsite2, ...]}
            for atom, vsite in vsite_pairs.items():
                _exclusions_vsite[vsite] = [atom]
                drude = drude_pairs.get(atom)
                if drude is not None:
                    _exclusions_vsite[vsite].append(drude)
            for atom1, atom2 in pairs12 + pairs13:
                vsite1 = vsite_pairs.get(atom1)
                vsite2 = vsite_pairs.get(atom2)
                drude1 = drude_pairs.get(atom1)
                drude2 = drude_pairs.get(atom2)
                if vsite1 is not None:
                    _exclusions_vsite[vsite1].append(atom2)
                    if drude2 is not None:
                        _exclusions_vsite[vsite1].append(drude2)
                if vsite2 is not None:
                    _exclusions_vsite[vsite2].append(atom1)
                    if drude1 is not None:
                        _exclusions_vsite[vsite2].append(drude1)
                if None not in [vsite1, vsite2]:
                    _exclusions_vsite[vsite1].append(vsite2)

            _exclusions = {**_exclusions_drude, **_exclusions_vsite}
            if len(_exclusions) > 0:
                string += '\n[ exclusions ]\n'
                string += ';  ai    aj    ...\n'
                for a1, a2_list in _exclusions.items():
                    string += '%5d' % (a1.id_in_mol + 1)
                    for a2 in a2_list:
                        string += ' %5d' % (a2.id_in_mol + 1)
                    string += '\n'

            string += '\n[ angles ]\n'
            for angle in mol.angles:
                if angle in system.constrain_angles:
                    continue
                aterm = system.angle_terms[angle]
                a1, a2, a3 = angle.atom1, angle.atom2, angle.atom3
                if aterm.__class__ in (HarmonicAngleTerm, SDKAngleTerm):
                    string += '%6i %6i %6i %6i %12.6f %12.4f\n' % (
                        a1.id_in_mol + 1, a2.id_in_mol + 1, a3.id_in_mol + 1,
                        1, aterm.theta * RAD2DEG, aterm.k * 2)
                elif aterm.__class__ is QuarticAngleTerm:
                    string += '%6i %6i %6i %6i %12.6f %12.4f\n' % (
                        a1.id_in_mol + 1, a2.id_in_mol + 1, a3.id_in_mol + 1,
                        1, aterm.theta * RAD2DEG, aterm.k2 * 2)
                elif aterm.__class__ is HarmonicCosineAngleTerm:
                    string += '%6i %6i %6i %6i %12.6f %12.4f\n' % (
                        a1.id_in_mol + 1, a2.id_in_mol + 1, a3.id_in_mol + 1,
                        2, aterm.theta * RAD2DEG, aterm.k * 2 / math.sin(aterm.theta) ** 2)
                elif aterm.__class__ is LinearAngleTerm:
                    bond12, bond23 = angle.bonds
                    bterm12 = system.bond_terms[bond12]
                    bterm23 = system.bond_terms[bond23]
                    a, k_linear = aterm.calc_a_k_linear(bterm12.length, bterm23.length)
                    string += '%6i %6i %6i %6i %12.6f %12.4f\n' % (
                        a1.id_in_mol + 1, a2.id_in_mol + 1, a3.id_in_mol + 1,
                        9, a, k_linear * 2)
                else:
                    raise Exception('Unsupported angle term')

            string += '\n[ dihedrals ]\n'
            for dihedral in mol.dihedrals:
                dterm = system.dihedral_terms[dihedral]
                if dterm.__class__ == OplsDihedralTerm:
                    dterm = dterm.to_periodic_term()
                a1, a2, a3, a4 = dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4
                if dterm.__class__ == PeriodicDihedralTerm:
                    for phase in dterm.phases:
                        string += '%6i %6i %6i %6i %6i %8.2f %12.4f %4i\n' % (
                            a1.id_in_mol + 1, a2.id_in_mol + 1,
                            a3.id_in_mol + 1, a4.id_in_mol + 1,
                            9, phase.phi * RAD2DEG, phase.k, phase.n)
                    if len(dterm.phases) == 0:
                        string += '%6i %6i %6i %6i %6i %8.2f %12.4f %4i  ; no dihedral parameters\n' % (
                            a1.id_in_mol + 1, a2.id_in_mol + 1,
                            a3.id_in_mol + 1, a4.id_in_mol + 1,
                            9, 0.0, 0.0, 1)
                else:
                    raise Exception('Unsupported dihedral term')

            string += '\n[ dihedrals ]\n'
            for improper in mol.impropers:
                iterm = system.improper_terms[improper]
                a1, a2, a3, a4 = improper.atom1, improper.atom2, improper.atom3, improper.atom4
                if iterm.__class__ == OplsImproperTerm:
                    # be careful about the sequence of atoms in OPLS improper definition
                    string += '%6i %6i %6i %6i %6i %8.2f %12.4f %4i\n' % (
                        a2.id_in_mol + 1, a3.id_in_mol + 1,
                        a1.id_in_mol + 1, a4.id_in_mol + 1,
                        4, 180, iterm.k, 2)
                elif iterm.__class__ == HarmonicImproperTerm:
                    string += '%6i %6i %6i %6i %6i %8.2f %12.4f\n' % (
                        a1.id_in_mol + 1, a2.id_in_mol + 1,
                        a3.id_in_mol + 1, a4.id_in_mol + 1,
                        2, iterm.phi * RAD2DEG, iterm.k * 2)
                else:
                    raise Exception('Unsupported improper term')

        string += '\n[ system ]\n'

        string += '\n[ molecules ]\n'
        for i, (mol, count) in enumerate(mols_unique.items()):
            string += '%s_%i %6i\n' % (mol.name, i + 1, count)

        with open(top_out, 'wb') as f:
            f.write(string.encode())

    @staticmethod
    def _export_mdp(system: System, mdp_out='grompp.mdp'):
        r_cut = system.ff.vdw_cutoff
        tau_t = 0.2 if DrudePolarTerm in system.ff.polar_term_classes else 1.0
        string = f'''; Created by mstk
integrator      = sd
dt              = 0.002 ; ps
nsteps          = 1000000

nstxout         = 0
nstvout         = 0
nstfout         = 0
nstxout-compressed = 1000
compressed-x-grps  = System

cutoff-scheme   = Verlet
pbc             = xyz
; rlist           = {r_cut}
coulombtype     = PME
rcoulomb        = {r_cut}
vdwtype         = Cut-off
rvdw            = {r_cut}
DispCorr        = EnerPres

tcoupl          = no; v-rescale
tc_grps         = System
tau_t           = {tau_t}
ref_t           = 300

pcoupl          = berendsen ; parrinello-rahman
pcoupltype      = isotropic
tau_p           = 0.5; 5
compressibility = 4.5e-5
ref_p           = 1.0

gen_vel         = yes
gen_temp        = 300

constraints     = h-bonds
constraint-algorithm = LINCS
'''
        with open(mdp_out, 'wb') as f:
            f.write(string.encode())
