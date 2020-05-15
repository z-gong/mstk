import itertools

import numpy as np
from .system import System
from ..forcefield import ForceField
from ..forcefield.ffterm import *
from ..topology import Topology, Atom, UnitCell, Psf, Bond, Angle, Dihedral, Improper
from ..trajectory import Frame, Trajectory, Gro
from .. import logger


class LammpsExporter():
    def __init__(self):
        pass

    @staticmethod
    def export(system: System, data_out='data.lmp', in_out='in.lmp'):
        supported_terms = {LJ126Term,
                           HarmonicBondTerm,
                           HarmonicAngleTerm,
                           PeriodicDihedralTerm,
                           OplsImproperTerm,
                           DrudeTerm}
        unsupported = system.ff_classes - supported_terms
        if unsupported != set():
            raise Exception('Unsupported FF terms: %s'
                            % (', '.join(map(lambda x: x.__name__, unsupported))))

        top = system.topology
        ff = system.ff
        bond_types = list(ff.bond_terms.values())
        angle_types = list(ff.angle_terms.values())
        dihedral_types = list(ff.dihedral_terms.values())
        improper_types = list(ff.improper_terms.values())

        ### assign lammps types ####################################
        atype_drude: AtomType = None
        lmp_types_real: {str: Atom} = {}  # {'type': Atom}
        lmp_types_drude: {str: Atom} = {}  # {'type': Atom}
        lmp_types_parent: {str: Atom} = {}  # {'type': Atom}
        lmp_types_atype: {str: AtomType} = {}  # {'type': AtomType}

        for atom in top.atoms:
            if atom.is_drude:
                continue
            if atom.type not in lmp_types_real:
                lmp_types_real[atom.type] = atom
                lmp_types_atype[atom.type] = ff.atom_types[atom.type]
        drude_pairs: {Atom: Atom} = dict(top.get_drude_pairs())
        for parent, drude in drude_pairs.items():
            if atype_drude is None:
                atype_drude = system.ff.atom_types[drude.type]
            typ = 'DP_' + parent.type
            if typ not in lmp_types_drude:
                lmp_types_drude[typ] = drude
                lmp_types_atype[typ] = atype_drude
            if parent.type not in lmp_types_parent:
                lmp_types_parent[parent.type] = parent

        lmp_type_list = list(lmp_types_real.keys()) + list(lmp_types_drude.keys())
        lmp_symbol_list = [atom.symbol or 'UNK' for atom in lmp_types_real.values()]
        fix_drude_list = ['C' if t in lmp_types_parent else 'D' if t in lmp_types_drude else 'N'
                          for t in lmp_type_list]
        lmp_shake_bonds = [i for i, bterm in enumerate(bond_types) if bterm.fixed]
        lmp_shake_angles = [i for i, aterm in enumerate(angle_types) if aterm.fixed]
        #############################################################

        string = 'Created by mstools\n'

        string += '\n%i atoms\n' % top.n_atom
        string += '%i bonds\n' % (top.n_bond + len(drude_pairs))
        string += '%i angles\n' % top.n_angle
        string += '%i dihedrals\n' % top.n_dihedral
        string += '%i impropers\n' % top.n_improper

        string += '\n%i atom types\n' % len(lmp_types_real)
        string += '%i bond types\n' % (len(bond_types) + len(lmp_types_parent))
        string += '%i angle types\n' % len(angle_types)
        string += '%i dihedral types\n' % len(dihedral_types)
        string += '%i improper types\n' % len(improper_types)

        if not top.cell.is_rectangular:
            raise Exception('Triclinic box haven\'t been implemented')
        box = top.cell.size
        string += '\n0 %10.4f  xlo xhi\n' % box[0]
        string += '0 %10.4f  ylo yhi\n' % box[1]
        string += '0 %10.4f  zlo zhi\n' % box[2]

        string += '\nMasses\n\n'

        for i, (typ, atom) in enumerate(lmp_types_real.items()):
            string += '%4i %8.3f  # %8s\n' % (i + 1, atom.mass, typ)
        for i, (typ, atom) in enumerate(lmp_types_drude.items()):
            string += '%4i %8.3f  # %8s\n' % (i + len(lmp_types_real) + 1, atom.mass, typ)

        string += '\nBond Coeffs  # harmonic\n\n'
        for i, bterm in enumerate(bond_types):
            string += '%4i %12.6f %10.4f  # %s-%s\n' % (
                i + 1, bterm.k / 4.184 / 100, bterm.length * 10, bterm.type1, bterm.type2)
        for i, parent in enumerate(lmp_types_parent.values()):
            pterm = system._polarizable_terms[parent]
            string += '%4i %12.6f %10.4f  # %s-%s\n' % (
                i + 1 + len(bond_types), pterm.k / 4.184 / 100, 0, parent.type, 'DP')

        string += '\nAngle Coeffs  # harmonic\n\n'
        for i, atype in enumerate(angle_types):
            string += '%4i %12.6f %10.4f  # %s-%s-%s\n' % (
                i + 1, atype.k / 4.184, atype.theta, atype.type1, atype.type2, atype.type3)

        string += '\nDihedral Coeffs  # opls\n\n'
        for i, dterm in enumerate(dihedral_types):
            dterm: PeriodicDihedralTerm
            if not dterm.follow_opls_convention:
                raise Exception('Only OPLS type dihedral is implemented')
            k1, k2, k3, k4 = map(lambda x: x * 2 / 4.184, dterm.get_opls_parameters())
            string += '%4i %12.6f %12.6f %12.6f %12.6f  # %s-%s-%s-%s\n' % (
                i + 1, k1, k2, k3, k4,
                dterm.type1, dterm.type2, dterm.type3, dterm.type4)

        string += '\nImproper Coeffs  # cvff\n\n'
        for i, iterm in enumerate(improper_types):
            string += '%4i %12.6f %4i %4i  # %s-%s-%s-%s\n' % (
                i + 1, iterm.k / 4.184, -1, 2, iterm.type2, iterm.type3, iterm.type1, iterm.type4)

        string += '\nAtoms\n\n'

        for atom in top.atoms:
            if atom.is_drude:
                lmp_type = 'DP_' + atom.bond_partners[0].type
            else:
                lmp_type = atom.type
            x, y, z = atom.position * 10
            string += '%8i %6i %6i %12.6f %10.4f %10.4f %10.4f  # %8s %8s\n' % (
                atom.id + 1, atom.molecule.id + 1, lmp_type_list.index(lmp_type) + 1,
                atom.charge, x, y, z, atom.name, atom.molecule.name)

        string += '\nBonds\n\n'

        for i, bond in enumerate(top.bonds):
            if not bond.is_drude:
                btype = bond_types.index(system._bond_terms[id(bond)]) + 1
            else:
                parent = bond.atom1 if bond.atom2.is_drude else bond.atom2
                btype = list(lmp_types_parent.keys()).index(parent.type) + len(bond_types) + 1
            string += '%6i %6i %6i %6i  # %s-%s\n' % (
                i + 1, btype, bond.atom1.id + 1, bond.atom2.id + 1, bond.atom1.name,
                bond.atom2.name)

        string += '\nAngles\n\n'

        for i, angle in enumerate(top.angles):
            atype = angle_types.index(system._angle_terms[id(angle)]) + 1
            atype1, atype2, a3 = angle.atom1, angle.atom2, angle.atom3
            string += '%6i %6i %6i %6i %6i  # %s-%s-%s\n' % (
                i + 1, atype, atype1.id + 1, atype2.id + 1, a3.id + 1, atype1.name, atype2.name,
                a3.name)

        string += '\nDihedrals\n\n'

        for i, dihedral in enumerate(top.dihedrals):
            dtype = dihedral_types.index(system._dihedral_terms[id(dihedral)]) + 1
            atype1, atype2, a3, a4 = dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4
            string += '%6i %6i %6i %6i %6i %6i  # %s-%s-%s-%s\n' % (
                i + 1, dtype, atype1.id + 1, atype2.id + 1, a3.id + 1, a4.id + 1,
                atype1.name, atype2.name, a3.name, a4.name)

        string += '\nImpropers\n\n'

        for i, improper in enumerate(top.impropers):
            itype = improper_types.index(system._improper_terms[id(improper)]) + 1
            atype1, atype2, a3, a4 = improper.atom1, improper.atom2, improper.atom3, improper.atom4
            string += '%6i %6i %6i %6i %6i %6i  # %s-%s-%s-%s\n' % (
                i + 1, itype, atype2.id + 1, a3.id + 1, atype1.id + 1, a4.id + 1,
                atype2.name, a3.name, atype1.name, a4.name)

        with open(data_out, 'w')  as f:
            f.write(string)

        lj14 = system.ff.scale_14_vdw
        coul14 = system.ff.scale_14_coulomb
        mix = 'geometric' if system.ff.lj_mixing_rule == ForceField.LJ_MIXING_GEOMETRIC else 'arthimatic'

        string = f'''# created by mstools
units real
boundary p p p
atom_style full
bond_style harmonic
angle_style harmonic
dihedral_style opls
improper_style cvff
special_bonds lj 0 0 {lj14} coul 0 0 {coul14}
'''

        cmd_shake = ''
        if len(lmp_shake_bonds) > 0:
            cmd_shake = 'fix SHAKE all shake b ' + ' '.join(map(str, lmp_shake_bonds))
            if len(lmp_shake_angles) > 0:
                cmd_shake += 'a ' + ' '.join(map(str, lmp_shake_angles))

        if DrudeTerm not in system.ff_classes:
            cmd_pair = ''
            for i, (typ_i, atom) in enumerate(lmp_types_real.items()):
                for j, (typ_j, atom) in enumerate(lmp_types_real.items()):
                    if j < i:
                        continue
                    atype1 = lmp_types_atype[typ_i]
                    atype2 = lmp_types_atype[typ_j]
                    try:
                        vdw = ff.get_vdw_term(atype1, atype2, mixing=False)
                    except:
                        continue
                    cmd_pair += 'pair_coeff %3i %3i %9.5f %8.4f  # %8s %8s %s\n' % (
                        i + 1, j + 1, vdw.epsilon / 4.184, vdw.sigma * 10,
                        lmp_type_list[i], lmp_type_list[j], ','.join(vdw.comments))

            string += f'''
pair_style lj/cut/coul/long 12.0
pair_modify mix {mix} tail yes
kspace_style pppm 1.0e-4

read_data data.lmp

{cmd_pair}

variable T equal 300
variable P equal 1
variable elec equal ecoul+elong

# thermo_style custom step temp press pe evdwl v_elec emol
# thermo 10
# minimize 1.0e-4 1.0e-6 1000 1000
# reset_timestep 0

{cmd_shake}
fix ICECUBE all momentum 100 linear 1 1 1

velocity all create $T 12345

timestep 1.0

fix NPT all npt temp $T $T 100 iso $P $P 1000

thermo_style custom step cpu temp press pe evdwl v_elec density
thermo 100
'''
        else:
            cmd_pair = ''
            for i, (typ_i, atom) in enumerate(lmp_types_real.items()):
                for j, (typ_j, atom) in enumerate(lmp_types_real.items()):
                    if j < i:
                        continue
                    atype1 = lmp_types_atype[typ_i]
                    atype2 = lmp_types_atype[typ_j]
                    try:
                        vdw = ff.get_vdw_term(atype1, atype2, mixing=False)
                    except:
                        continue
                    cmd_pair += 'pair_coeff %3i %3i %9.5f %8.4f %8.4f %5.3f  # %8s %8s %s\n' % (
                        i + 1, j + 1, vdw.epsilon / 4.184, vdw.sigma * 10,
                        atom._alpha * 1000, atom._thole,
                        lmp_type_list[i], lmp_type_list[j], ','.join(vdw.comments))
            for i, (typ, atom) in enumerate(lmp_types_parent.items()):
                ii = i + len(lmp_types_real)
                vdw = ff.get_vdw_term(atype_drude, atype_drude)
                cmd_pair += 'pair_coeff %3i %3i %9.5f %8.4f %8.4f %5.3f  # %8s %8s %s\n' % (
                    ii + 1, ii + 1, vdw.epsilon / 4.184, vdw.sigma * 10,
                    atom._alpha * 1000, atom._thole,
                    lmp_type_list[ii], lmp_type_list[ii], ','.join(vdw.comments))

            string += f'''
pair_style lj/cut/thole/long 2.6 12.0
pair_modify mix {mix} tail yes
kspace_style pppm 1.0e-4

read_data data.lmp extra/special/per/atom 99

{cmd_pair}

group ATOMS type {' '.join(map(str, [i + 1 for i, x in enumerate(fix_drude_list) if x != 'D']))}
group CORES type {' '.join(map(str, [i + 1 for i, x in enumerate(fix_drude_list) if x == 'C']))}
group DRUDES type {' '.join(map(str, [i + 1 for i, x in enumerate(fix_drude_list) if x == 'D']))}
fix DRUDE all drude {' '.join(fix_drude_list)}

variable T equal 300
variable P equal 1
variable elec equal e_coul+e_long

# thermo_style custom step temp press pe evdwl v_elec emol
# thermo 10
# minimize 1.0e-4 1.0e-6 1000 1000
# reset_timestep 0

{cmd_shake}
fix ICECUBE all momentum 100 linear 1 1 1

velocity all create $T 12345

comm_modify vel yes
compute TDRUDE all temp/drude

timestep 1.0

fix SD  all langevin/drude $TK 200 12345 1 50 23456
fix NPH all nph iso $P $P 1000
# fix NPT all tgnpt/drude temp $T $T 100 1 25 iso $P $P 1000

thermo_style custom step cpu c_TDRUDE[3] c_TDRUDE[4] c_TDRUDE[2] press pe evdwl v_elec density
thermo 100
'''

        string += f'''
dump TRAJ all custom 10000 dump.lammpstrj id mol type element q xu yu zu
dump_modify TRAJ sort id element {' '.join(lmp_symbol_list)}
dump DCD all dcd 1000 dump.dcd
dump_modify DCD unwrap yes

run 1000000
write_restart rst
'''

        with open(in_out, 'w') as f:
            f.write(string)
