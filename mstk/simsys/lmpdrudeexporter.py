import os
from .system import System
from ..forcefield import *
from ..topology import *
from .. import logger


class LammpsDrudeExporter():
    '''
    LammpsDrudeExporter export a Drude polarizable :class:`System` to input files for LAMMPS.

    The following potential functions are currently supported:

    * :class:`~mstk.forcefield.LJ126Term`
    * :class:`~mstk.forcefield.HarmonicBondTerm`
    * :class:`~mstk.forcefield.HarmonicAngleTerm`
    * :class:`~mstk.forcefield.PeriodicDihedralTerm`
    * :class:`~mstk.forcefield.OplsImproperTerm`
    * :class:`~mstk.forcefield.DrudeTerm`

    In order to run simulation of Drude model, the package DRUDE should be included when compiling Lammps.
    '''

    def __init__(self):
        pass

    @staticmethod
    def export(system, data_out, in_out, **kwargs):
        '''
        Generate input files for Lammps from a system

        Parameters
        ----------
        system : System
        data_out : str
        in_out : str
        '''
        if not system.use_pbc:
            raise Exception('PBC required for exporting LAMMPS')

        supported_terms = {LJ126Term,
                           HarmonicBondTerm,
                           HarmonicAngleTerm, LinearAngleTerm,
                           PeriodicDihedralTerm,
                           OplsImproperTerm,
                           DrudeTerm}
        unsupported = system.ff_classes - supported_terms
        if unsupported != set():
            raise Exception('Unsupported FF terms: %s' % (', '.join(map(lambda x: x.__name__, unsupported))))

        if system.topology.has_virtual_site:
            raise Exception('Virtual sites not supported by LAMMPS')

        if LinearAngleTerm in system.ff_classes:
            logger.warning('LinearAngleTerm not supported by LAMMPS. Exported in harmonic form at 178 degree')

        top = system.topology
        ff = system.ff

        ### Assign LAMMPS types ####################################
        bond_types = list(ff.bond_terms.values())  # bond types between Drude pairs not included
        angle_types = list(ff.angle_terms.values())
        dihedral_types = list(ff.dihedral_terms.values())
        improper_types = list(ff.improper_terms.values())

        # In LAMMPS, Drude parameters are assigned with pair coefficients
        # Drude particles attached to different heavy atoms must have its own LAMMPS atom type
        lmp_types_real: {str: Atom} = {}  # {'typ': Atom}
        lmp_types_drude: {str: Atom} = {}  # {'typ': Atom}
        lmp_types_parent: {str: Atom} = {}  # {'typ': Atom}
        lmp_types_atype: {str: AtomType} = {}  # {'typ': AtomType}
        atype_drude: AtomType = None

        for atom in top.atoms:
            if atom.is_drude:
                continue
            if atom.type not in lmp_types_real:
                lmp_types_real[atom.type] = atom
                lmp_types_atype[atom.type] = ff.atom_types[atom.type]
        for parent, drude in system.drude_pairs.items():
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
        lmp_symbol_list += [atom.symbol or 'UNK' for atom in lmp_types_drude.values()]
        fix_drude_list = ['C' if t in lmp_types_parent else 'D' if t in lmp_types_drude else 'N'
                          for t in lmp_type_list]
        lmp_shake_bonds = [i + 1 for i, bterm in enumerate(bond_types) if bterm.fixed]
        lmp_shake_angles = [i + 1 for i, aterm in enumerate(angle_types) if aterm.fixed]
        #############################################################

        string = 'Created by mstk\n'

        string += '\n%i atoms\n' % top.n_atom
        string += '%i bonds\n' % top.n_bond
        string += '%i angles\n' % top.n_angle
        string += '%i dihedrals\n' % top.n_dihedral
        string += '%i impropers\n' % top.n_improper

        string += '\n%i atom types\n' % len(lmp_type_list)
        string += '%i bond types\n' % (len(bond_types) + len(lmp_types_parent))
        string += '%i angle types\n' % len(angle_types)
        string += '%i dihedral types\n' % len(dihedral_types)
        string += '%i improper types\n' % len(improper_types)

        if not top.cell.is_rectangular:
            raise Exception('Triclinic box haven\'t been implemented')
        box = top.cell.size * 10
        string += '\n0 %10.4f  xlo xhi\n' % box[0]
        string += '0 %10.4f  ylo yhi\n' % box[1]
        string += '0 %10.4f  zlo zhi\n' % box[2]

        string += '\nMasses\n\n'

        for i, (typ, atom) in enumerate(lmp_types_real.items()):
            string += '%4i %8.3f  # %8s\n' % (i + 1, atom.mass, typ)
        for i, (typ, atom) in enumerate(lmp_types_drude.items()):
            string += '%4i %8.3f  # %8s\n' % (i + len(lmp_types_real) + 1, atom.mass, typ)

        if len(bond_types) + len(lmp_types_parent) > 0:
            string += '\nBond Coeffs  # harmonic\n\n'

        for i, bterm in enumerate(bond_types):
            string += '%4i %12.6f %10.4f  # %s-%s\n' % (
                i + 1, bterm.k / 4.184 / 100, bterm.length * 10, bterm.type1, bterm.type2)
        for i, parent in enumerate(lmp_types_parent.values()):
            pterm = system.polarizable_terms[parent]
            string += '%4i %12.6f %10.4f  # %s-%s\n' % (
                i + 1 + len(bond_types), pterm.k / 4.184 / 100, 0, parent.type, 'DP')

        if len(angle_types) > 0:
            string += '\nAngle Coeffs  # harmonic\n\n'

        for i, atype in enumerate(angle_types):
            string += '%4i %12.6f %10.4f  # %s-%s-%s\n' % (
                i + 1, atype.k / 4.184, min(atype.theta, 178), atype.type1, atype.type2, atype.type3)

        if len(dihedral_types) > 0:
            string += '\nDihedral Coeffs  # opls\n\n'

        for i, dterm in enumerate(dihedral_types):
            dterm: PeriodicDihedralTerm
            if not dterm.is_opls_convention:
                raise Exception('Only OPLS type dihedral is implemented')
            k1, k2, k3, k4 = map(lambda x: x * 2 / 4.184, dterm.get_opls_parameters())
            string += '%4i %12.6f %12.6f %12.6f %12.6f  # %s-%s-%s-%s\n' % (
                i + 1, k1, k2, k3, k4,
                dterm.type1, dterm.type2, dterm.type3, dterm.type4)

        if len(improper_types) > 0:
            string += '\nImproper Coeffs  # cvff\n\n'

        for i, iterm in enumerate(improper_types):
            string += '%4i %12.6f %4i %4i  # %s-%s-%s-%s\n' % (
                i + 1, iterm.k / 4.184, -1, 2, iterm.type2, iterm.type3, iterm.type1, iterm.type4)

        string += '\nAtoms\n\n'

        for atom in top.atoms:
            if atom.is_drude:
                typ = 'DP_' + atom.bond_partners[0].type
            else:
                typ = atom.type
            x, y, z = atom.position * 10
            string += '%8i %6i %6i %12.6f %10.4f %10.4f %10.4f  # %8s %8s\n' % (
                atom.id + 1, atom.molecule.id + 1, lmp_type_list.index(typ) + 1,
                atom.charge, x, y, z, atom.name, atom.molecule.name)

        if top.n_bond > 0:
            string += '\nBonds\n\n'

        for i, bond in enumerate(top.bonds):
            a1, a2 = bond.atom1, bond.atom2
            if not bond.is_drude:
                btype = bond_types.index(system.bond_terms[id(bond)]) + 1
            else:
                parent = a1 if a2.is_drude else a2
                btype = list(lmp_types_parent.keys()).index(parent.type) + len(bond_types) + 1
            string += '%6i %6i %6i %6i  # %s-%s\n' % (
                i + 1, btype, a1.id + 1, a2.id + 1, a1.name, a2.name)

        if top.n_angle > 0:
            string += '\nAngles\n\n'

        for i, angle in enumerate(top.angles):
            atype = angle_types.index(system.angle_terms[id(angle)]) + 1
            a1, a2, a3 = angle.atom1, angle.atom2, angle.atom3
            string += '%6i %6i %6i %6i %6i  # %s-%s-%s\n' % (
                i + 1, atype, a1.id + 1, a2.id + 1, a3.id + 1, a1.name, a2.name, a3.name)

        if top.n_dihedral > 0:
            string += '\nDihedrals\n\n'

        for i, dihedral in enumerate(top.dihedrals):
            dtype = dihedral_types.index(system.dihedral_terms[id(dihedral)]) + 1
            a1, a2, a3, a4 = dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4
            string += '%6i %6i %6i %6i %6i %6i  # %s-%s-%s-%s\n' % (
                i + 1, dtype, a1.id + 1, a2.id + 1, a3.id + 1, a4.id + 1,
                a1.name, a2.name, a3.name, a4.name)

        if top.n_improper > 0:
            string += '\nImpropers\n\n'

        for i, improper in enumerate(top.impropers):
            itype = improper_types.index(system.improper_terms[id(improper)]) + 1
            a1, a2, a3, a4 = improper.atom1, improper.atom2, improper.atom3, improper.atom4
            string += '%6i %6i %6i %6i %6i %6i  # %s-%s-%s-%s\n' % (
                i + 1, itype, a2.id + 1, a3.id + 1, a1.id + 1, a4.id + 1,
                a2.name, a3.name, a1.name, a4.name)

        with open(data_out, 'wb')  as f:
            f.write(string.encode())

        cmd_mix = 'geometric'
        if system.ff.lj_mixing_rule == ForceField.LJ_MIXING_LB:
            cmd_mix = 'arithmetic'

        string = f'''# created by mstk
units real
boundary p p p
atom_style full
bond_style harmonic
angle_style harmonic
dihedral_style opls
improper_style cvff
special_bonds lj 0 0 {ff.scale_14_vdw} coul 0 0 {ff.scale_14_coulomb}
'''

        cmd_pair = ''
        for i, (typ_i, a1) in enumerate(lmp_types_real.items()):
            for j, (typ_j, a2) in enumerate(lmp_types_real.items()):
                if j < i:
                    continue
                atype1 = lmp_types_atype[typ_i]
                atype2 = lmp_types_atype[typ_j]
                try:
                    vdw = ff.get_vdw_term(atype1, atype2, mixing=False)
                except:
                    continue
                cmd_pair += 'pair_coeff %3i %3i %9.5f %8.4f %8.4f %7.3f  # %8s %8s %s\n' % (
                    i + 1, j + 1, vdw.epsilon / 4.184, vdw.sigma * 10,
                    math.sqrt(a1.alpha * a2.alpha) * 1000, (a1.thole + a2.thole) / 2,
                    lmp_type_list[i], lmp_type_list[j], ','.join(vdw.comments))
        for i, (typ, atom) in enumerate(lmp_types_parent.items()):
            ii = i + len(lmp_types_real)
            vdw = ff.get_vdw_term(atype_drude, atype_drude)
            cmd_pair += 'pair_coeff %3i %3i %9.5f %8.4f %8.4f %7.3f  # %8s %8s %s\n' % (
                ii + 1, ii + 1, vdw.epsilon / 4.184, vdw.sigma * 10,
                atom.alpha * 1000, atom.thole,
                lmp_type_list[ii], lmp_type_list[ii], ','.join(vdw.comments))

        string += f'''
pair_style lj/cut/thole/long 2.6 12.0
pair_modify mix {cmd_mix} tail yes
kspace_style pppm 1.0e-4

read_data {os.path.basename(data_out)} extra/special/per/atom 99

{cmd_pair}
group ATOMS type {' '.join(map(str, [i + 1 for i, x in enumerate(fix_drude_list) if x != 'D']))}
group CORES type {' '.join(map(str, [i + 1 for i, x in enumerate(fix_drude_list) if x == 'C']))}
group DRUDES type {' '.join(map(str, [i + 1 for i, x in enumerate(fix_drude_list) if x == 'D']))}
fix DRUDE all drude {' '.join(fix_drude_list)}
'''

        cmd_shake = ''
        if len(lmp_shake_bonds) > 0:
            cmd_shake = 'fix SHAKE all shake 0.0001 20 0 b ' + ' '.join(map(str, lmp_shake_bonds))
            if len(lmp_shake_angles) > 0:
                cmd_shake += ' a ' + ' '.join(map(str, lmp_shake_angles))
        string += f'''
variable T equal 300
variable P equal 1
variable elec equal ecoul+elong

# thermo_style custom step press pe ebond eangle edihed eimp evdwl v_elec
# thermo 10
# minimize 1.0e-4 1.0e-6 200 1000
# reset_timestep 0

{cmd_shake}
fix ICECUBE all momentum 100 linear 1 1 1

velocity all create $T 12345
timestep 1.0

comm_modify vel yes
compute TDRUDE all temp/drude

fix NPT all tgnpt/drude temp $T $T 100 1 25 iso $P $P 1000 fixedpoint 0 0 0
# fix SD  all langevin/drude $T 200 12345 1 50 23456
# fix NPH all nph iso $P $P 1000

thermo_style custom step cpu c_TDRUDE[3] c_TDRUDE[4] c_TDRUDE[2] press pe emol evdwl v_elec density
thermo_modify flush yes
thermo 1000

variable slog equal logfreq(10,9,10)
dump TRJ all custom 10 dump.lammpstrj id mol type element q xu yu zu
dump_modify TRJ sort id element {' '.join(lmp_symbol_list)} every v_slog first yes
dump DCD all dcd 10000 dump.dcd
dump_modify DCD unwrap yes

restart 1000000 rst_*
run 1000000
'''

        with open(in_out, 'wb') as f:
            f.write(string.encode())
