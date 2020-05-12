import itertools
import numpy as np
from .system import System
from ..forcefield.ffterm import *
from ..forcefield import FFSet
from ..topology import Topology, Atom, UnitCell, Psf, Bond, Angle, Dihedral, Improper, Molecule
from ..trajectory import Frame, Trajectory, Gro
from .. import logger


class GromacsExporter():
    def __init__(self):
        pass

    @staticmethod
    def export(system: System, gro_out, top_out, mdp_out):
        if gro_out is not None:
            GromacsExporter.export_gro(system, gro_out)
        if top_out is not None:
            GromacsExporter.export_top(system, top_out)
        if mdp_out is not None:
            GromacsExporter.export_mdp(mdp_out)

    @staticmethod
    def export_gro(system, gro_out='conf.gro'):
        frame = Frame(system._topology.n_atom)
        frame.cell = system._topology.cell
        frame.positions = system._topology.positions
        gro = Gro(gro_out, 'w')
        gro.write_frame(frame, system._topology)
        gro.close()

    @staticmethod
    def export_top(system, top_out='topol.top'):
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
                           'The exported files should only be used for analyzing trajectory')
        if SDKAngleTerm in system.ff_classes:
            logger.warning('SDKAngleTerm not supported by GROMACS, '
                           'replaced by harmonic function')
        if DrudeTerm in system.ff_classes:
            logger.warning('DurdeTerm havn\'t been implemented for GROMACS. '
                           'Will treat Drude particles as normal atoms')

        mols_unique = system._topology.get_unique_molecules()

        string = '; GROMACS topol file created by mstools\n'

        string += '\n[ defaults ]\n'
        string += '; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n'
        _comb = 2 if system.ff.lj_mixing_rule == FFSet.LJ_MIXING_LB else 3
        string += '%6i %6i %12s %10.4f %10.4f\n' % (
            1, _comb, 'yes', system.ff.scale_14_vdw, system.ff.scale_14_coulomb)

        string += '\n[ atomtypes ]\n'
        string += ';     name       mass     charge      ptype      sigma      epsilon\n'
        for atype in system.ff.atom_types.values():
            vdw = system.ff.get_vdw_term(atype, atype)
            if vdw.__class__ in (LJ126Term, MieTerm):
                string += '%10s %10.4f %12.6f %6s %12.6f %12.6f  ; %s\n' % (
                    atype.name, 0.0, 0.0, 'A', vdw.sigma, vdw.epsilon, ' '.join(vdw.comments))
            else:
                raise Exception('Unsupported vdW term')

        string += '\n[ nonbond_params ]\n'
        string += ';       i        j        func      sigma      epsilon\n'
        for at1, at2 in itertools.combinations(system.ff.atom_types.values(), 2):
            vdw = system.ff.get_vdw_term(at1, at2, mixing=False)
            if vdw is not None:
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
                string += '%6i %10s %6i %10s %6s %6i %12.6f %10.4f\n' % (
                    j + 1, atom.type, 1, mol.name, atom.symbol, 1, atom.charge, atom.mass)

            string += '\n[ pairs ]\n'
            _, _, pair14 = mol.get_12_13_14_pairs()
            for a1, a2 in pair14:
                string += '%6i %6i %6i\n' % (a1._id_in_molecule + 1, a2._id_in_molecule + 1, 1)

            string += '\n[ constraints ]\n'
            for bond in mol.bonds:
                distance = system._constrain_bonds.get(id(bond))
                if distance is not None:
                    a1, a2 = bond.atom1, bond.atom2
                    string += '%6i %6i %6i %12.6f\n' % (
                        a1._id_in_molecule + 1, a2._id_in_molecule + 1, 1, distance)
            for angle in mol.angles:
                distance = system._constrain_angles.get(id(angle))
                if distance is not None:
                    a1, a3 = angle.atom1, angle.atom3
                    string += '%6i %6i %6i %12.6f\n' % (
                        a1._id_in_molecule + 1, a3._id_in_molecule + 1, 2, distance)

            string += '\n[ bonds ]\n'
            for bond in mol.bonds:
                if bond.is_drude:
                    # TODO
                    continue
                bterm = system._bond_terms[id(bond)]
                if bterm.__class__ == HarmonicBondTerm:
                    a1, a2 = bond.atom1, bond.atom2
                    string += '%6i %6i %6i %12.6f %12.4f\n' % (
                        a1._id_in_molecule + 1, a2._id_in_molecule + 1,
                        1, bterm.length, bterm.k * 2)
                else:
                    raise Exception('Unsupported bond term')

            string += '\n[ angles ]\n'
            for angle in mol.angles:
                aterm = system._angle_terms[id(angle)]
                a1, a2, a3 = angle.atom1, angle.atom2, angle.atom3
                if aterm.__class__ in (HarmonicAngleTerm, SDKAngleTerm):
                    string += '%6i %6i %6i %6i %12.6f %12.4f\n' % (
                        a1._id_in_molecule + 1, a2._id_in_molecule + 1, a3._id_in_molecule + 1,
                        1, aterm.theta, aterm.k * 2)
                else:
                    raise Exception('Unsupported angle term')

            string += '\n[ dihedrals ]\n'
            for dihedral in mol.dihedrals:
                dterm = system._dihedral_terms[id(dihedral)]
                a1, a2, a3, a4 = dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4
                if dterm.__class__ == PeriodicDihedralTerm:
                    for para in dterm.parameters:
                        string += '%6i %6i %6i %6i %6i %.2f %12.4f %4i\n' % (
                            a1._id_in_molecule + 1, a2._id_in_molecule + 1,
                            a3._id_in_molecule + 1, a4._id_in_molecule + 1,
                            9, para.phi, para.k, para.n)
                else:
                    raise Exception('Unsupported dihedral term')

            string += '\n[ dihedrals ]\n'
            for improper in mol.impropers:
                iterm = system._improper_terms[id(improper)]
                a1, a2, a3, a4 = improper.atom1, improper.atom2, improper.atom3, improper.atom4
                if iterm.__class__ == OplsImproperTerm:
                    # be careful about the sequence of atoms in OPLS improper definition
                    string += '%6i %6i %6i %6i %6i %.2f %12.4f %4i\n' % (
                        a2._id_in_molecule + 1, a3._id_in_molecule + 1,
                        a1._id_in_molecule + 1, a4._id_in_molecule + 1,
                        4, 180, iterm.k, 2)
                elif iterm.__class__ == HarmonicImproperTerm:
                    string += '%6i %6i %6i %6i %6i %.2f %12.4f\n' % (
                        a1._id_in_molecule + 1, a2._id_in_molecule + 1,
                        a3._id_in_molecule + 1, a4._id_in_molecule + 1,
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
    def export_mdp(mdp_out='grompp.mdp'):
        with open(mdp_out, 'w')  as f:
            f.write('''; Created by mstools
integrator      = sd
dt              = 0.002 ; ps
nsteps          = 1000000

nstxout         = 0
nstvout         = 0
nstfout         = 0
nstxout-compressed = 1000
compressed-x-grps  = System

cutoff-scheme   = Verlet
rlist           = 1.2
coulombtype     = PME
rcoulomb        = 1.2
rvdw            = 1.2
DispCorr        = EnerPres

Tcoupl          = no; v-rescale
tc_grps         = System
tau_t           = 1.0
ref_t           = 300

Pcoupl          = berendsen ; parrinello-rahman
pcoupltype      = isotropic
tau_p           = 0.5; 5
compressibility = 4.5e-5
ref_p           = 1.0

gen_vel         = yes
gen_temp        = 300

constraints     = h-bonds
constraint-algorithm = LINCS
''')
