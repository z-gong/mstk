import itertools
from .system import System
from ..forcefield import *
from ..topology import *
from ..trajectory import *
from .. import logger


class GromacsExporter():
    '''
    GromacsExporter export a :class:`System` to input files for Gromacs.

    The following potential functions are currently supported:

    * :class:`~mstools.forcefield.LJ126Term`
    * :class:`~mstools.forcefield.HarmonicBondTerm`
    * :class:`~mstools.forcefield.HarmonicAngleTerm`
    * :class:`~mstools.forcefield.PeriodicDihedralTerm`
    * :class:`~mstools.forcefield.OplsImproperTerm`
    * :class:`~mstools.forcefield.HarmonicImproperTerm`
    * :class:`~mstools.forcefield.DrudeTerm`

    The :class:`~mstools.forcefield.MieTerm` can be exported, but it will be in the LJ-12-6 form.
    If the dispersion of `MieTerm` is in 6-th power, then MieTerm can be correctly handled
    by small modifications to topol and mdp files after exported.
    Refer to the GROMACS documentation for details.

    The :class:`~mstools.forcefield.SDKAngleTerm` can be exported, but it will be in the harmonic form.
    Usually it is acceptable if the molecule is not going to form a coil structure.
    Refer to the original SDK paper for details.
    '''

    def __init__(self):
        pass

    @staticmethod
    def export(system, gro_out, top_out, mdp_out):
        '''
        Generate input files for Gromacs from a system

        Parameters
        ----------
        system : System
        gro_out : str or None
        top_out : str or None
        mdp_out : str or None
        '''
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

        if MieTerm in system.ff_classes:
            logger.warning('MieTerm not supported by GROMACS. '
                           'Will be exported in LJ-126 form')
        if SDKAngleTerm in system.ff_classes:
            logger.warning('SDKAngleTerm not supported by GROMACS. '
                           'Will be exported in harmonic form')

        mols_unique = system._topology.get_unique_molecules(deepcopy=False)

        string = '; GROMACS topol file created by mstools\n'

        string += '\n[ defaults ]\n'
        string += '; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n'
        _comb = 2 if system.ff.lj_mixing_rule == ForceField.LJ_MIXING_LB else 3
        string += '%6i %6i %12s %10.4f %10.4f\n' % (
            1, _comb, 'yes', system.ff.scale_14_vdw, system.ff.scale_14_coulomb)

        string += '\n[ atomtypes ]\n'
        string += ';     name       mass     charge      ptype      sigma      epsilon\n'
        drude_types = set()
        if DrudeTerm in system.ff_classes:
            for atom in system.topology.atoms:
                if atom.is_drude:
                    drude_types.add(atom.type)
        for atype in system.ff.atom_types.values():
            ptype = 'S' if atype.name in drude_types else 'A'
            vdw = system.ff.get_vdw_term(atype, atype)
            if vdw.__class__ in (LJ126Term, MieTerm):
                string += '%10s %10.4f %12.6f %6s %12.6f %12.6f  ; %s\n' % (
                    atype.name, 0.0, 0.0, ptype, vdw.sigma, vdw.epsilon, ' '.join(vdw.comments))
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
                    j + 1, atom.type, 1, mol.name, atom.symbol, j + 1, atom.charge, mass)

            string += '\n[ pairs ]\n'
            pairs12, pairs13, pairs14 = mol.get_12_13_14_pairs()
            for a1, a2 in pairs14:
                string += '%6i %6i %6i\n' % (a1.id_in_mol + 1, a2.id_in_mol + 1, 1)

            if DrudeTerm in system.ff_classes:
                drude_pairs = dict(mol.get_drude_pairs())
                _pairs14_drude = []
                for a1, a2 in pairs14:
                    if a1 in drude_pairs:
                        _pairs14_drude.append((drude_pairs[a1], a2))
                    if a2 in drude_pairs:
                        _pairs14_drude.append((a1, drude_pairs[a2]))
                    if a1 in drude_pairs and a2 in drude_pairs:
                        _pairs14_drude.append((drude_pairs[a1], drude_pairs[a2]))
                for a1, a2 in _pairs14_drude:
                    string += '%6i %6i %6i\n' % (a1.id_in_mol + 1, a2.id_in_mol + 1, 1)

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
                distance = system.constrain_bonds.get(id(bond))
                if distance is not None:
                    a1, a2 = bond.atom1, bond.atom2
                    string += '%6i %6i %6i %12.6f\n' % (
                        a1.id_in_mol + 1, a2.id_in_mol + 1, 1, distance)
            for angle in mol.angles:
                distance = system.constrain_angles.get(id(angle))
                if distance is not None:
                    a1, a3 = angle.atom1, angle.atom3
                    string += '%6i %6i %6i %12.6f\n' % (
                        a1.id_in_mol + 1, a3.id_in_mol + 1, 2, distance)

            string += '\n[ bonds ]\n'
            for bond in mol.bonds:
                if bond.is_drude:
                    continue
                if id(bond) in system.constrain_bonds:
                    continue
                bterm = system.bond_terms[id(bond)]
                if bterm.__class__ == HarmonicBondTerm:
                    a1, a2 = bond.atom1, bond.atom2
                    string += '%6i %6i %6i %12.6f %12.4f\n' % (
                        a1.id_in_mol + 1, a2.id_in_mol + 1,
                        1, bterm.length, bterm.k * 2)
                else:
                    raise Exception('Unsupported bond term')

            if DrudeTerm in system.ff_classes:
                string += '\n[ exclusions ]\n'
                string += ';  ai    aj    ...\n'
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
                for d1 in _exclusions_drude.keys():
                    string += '%5d' % (d1.id_in_mol + 1)
                    for d2 in _exclusions_drude[d1]:
                        string += ' %5d' % (d2.id_in_mol + 1)
                    string += '\n'
                string += '\n'

            string += '\n[ angles ]\n'
            for angle in mol.angles:
                if id(angle) in system.constrain_angles:
                    continue
                aterm = system.angle_terms[id(angle)]
                a1, a2, a3 = angle.atom1, angle.atom2, angle.atom3
                if aterm.__class__ in (HarmonicAngleTerm, SDKAngleTerm):
                    string += '%6i %6i %6i %6i %12.6f %12.4f\n' % (
                        a1.id_in_mol + 1, a2.id_in_mol + 1, a3.id_in_mol + 1,
                        1, aterm.theta, aterm.k * 2)
                else:
                    raise Exception('Unsupported angle term')

            string += '\n[ dihedrals ]\n'
            for dihedral in mol.dihedrals:
                dterm = system.dihedral_terms[id(dihedral)]
                a1, a2, a3, a4 = dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4
                if dterm.__class__ == PeriodicDihedralTerm:
                    for para in dterm.parameters:
                        string += '%6i %6i %6i %6i %6i %8.2f %12.4f %4i\n' % (
                            a1.id_in_mol + 1, a2.id_in_mol + 1,
                            a3.id_in_mol + 1, a4.id_in_mol + 1,
                            9, para.phi, para.k, para.n)
                    if len(dterm.parameters) == 0:
                        string += '%6i %6i %6i %6i %6i %8.2f %12.4f %4i  ; no dihedral parameters\n' % (
                            a1.id_in_mol + 1, a2.id_in_mol + 1,
                            a3.id_in_mol + 1, a4.id_in_mol + 1,
                            9, 0.0, 0.0, 1)
                else:
                    raise Exception('Unsupported dihedral term')

            string += '\n[ dihedrals ]\n'
            for improper in mol.impropers:
                iterm = system.improper_terms[id(improper)]
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
                        2, iterm.phi, iterm.k * 2)
                else:
                    raise Exception('Unsupported improper term')

        string += '\n[ system ]\n'

        string += '\n[ molecules ]\n'
        for i, (mol, count) in enumerate(mols_unique.items()):
            string += '%s_%i %6i\n' % (mol.name, i + 1, count)

        with open(top_out, 'w') as f:
            f.write(string)

    @staticmethod
    def _export_mdp(system: System, mdp_out='grompp.mdp'):
        tau_t = 0.2 if DrudeTerm in system.ff_classes else 1.0
        with open(mdp_out, 'w')  as f:
            f.write(f'''; Created by mstools
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
; rlist           = 1.2
coulombtype     = PME
rcoulomb        = 1.2
vdwtype         = Cut-off
rvdw            = 1.2
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
''')
