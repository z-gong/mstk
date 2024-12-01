import os
from .system import System
from ..forcefield import *
from ..topology import *
from .. import logger


class LammpsExporter:
    '''
    LammpsExporter export a non-polarizable :class:`System` to input files for LAMMPS.

    LAMMPS is powerful and flexible. But the input file for LAMMPS is a mess and cannot be handled in an elegant way,
    especially when there are hybrid functional forms.
    Here, only limited types of FF (OPLS-AA and SDK-CG) will be considered.

    The following potential functions are currently supported:

    * :class:`~mstk.forcefield.LJ126Term`
    * :class:`~mstk.forcefield.MieTerm`
    * :class:`~mstk.forcefield.HarmonicBondTerm`
    * :class:`~mstk.forcefield.HarmonicAngleTerm`
    * :class:`~mstk.forcefield.SDKAngleTerm`
    * :class:`~mstk.forcefield.OplsDihedralTerm`
    * :class:`~mstk.forcefield.PeriodicDihedralTerm`
    * :class:`~mstk.forcefield.OplsImproperTerm`

    The :class:~mstk.forcefield.LinearAngleTerm can be exported, but it will be in the harmonic form.
    This is only valid when the force constant is large enough to keep the angle nearby 180 degree.

    In order to run simulation of SDK-CG model, the package CG-SDK should be included when compiling LAMMPS.
    '''

    def __init__(self):
        pass

    @staticmethod
    def export(system, data_out, in_out, **kwargs):
        '''
        Generate input files for LAMMPS from a system

        Parameters
        ----------
        system : System
        data_out : str
        in_out : str
        '''
        if not system.use_pbc:
            raise Exception('PBC required for exporting LAMMPS')

        supported_terms = {LJ126Term, MieTerm,
                           HarmonicBondTerm, MorseBondTerm,
                           HarmonicAngleTerm, SDKAngleTerm, LinearAngleTerm,
                           OplsDihedralTerm, PeriodicDihedralTerm,
                           OplsImproperTerm}
        unsupported = system.ff.energy_term_classes - supported_terms
        if unsupported != set():
            raise Exception('Unsupported FF terms: %s' % (', '.join(map(lambda x: x.__name__, unsupported))))

        if system.topology.has_virtual_site:
            raise Exception('Virtual sites not supported by LAMMPS')

        if LinearAngleTerm in system.ff.angle_term_classes:
            logger.warning('LinearAngleTerm not supported by LAMMPS. Exported in harmonic form')

        top = system.topology
        ff = system.ff

        # Morse bond requires hybrid bond styles
        is_morse = MorseBondTerm in system.ff.bond_term_classes

        # SDK CGFF requires all pairwise parameters and hybrid angle styles
        is_sdk = MieTerm in system.ff.vdw_term_classes or SDKAngleTerm in system.ff.angle_term_classes

        ### Assign LAMMPS types ####################################
        bond_types = list(ff.bond_terms.values())
        angle_types = list(ff.angle_terms.values())
        dihedral_types = list(ff.dihedral_terms.values())
        improper_types = list(ff.improper_terms.values())

        lmp_types: {str: Atom} = {}  # {'typ': Atom}
        lmp_types_atype: {str: AtomType} = {}  # {'typ': AtomType}

        for atom in top.atoms:
            if atom.type not in lmp_types:
                lmp_types[atom.type] = atom
                lmp_types_atype[atom.type] = ff.atom_types[atom.type]

        lmp_type_list = list(lmp_types.keys())
        lmp_symbol_list = [atom.symbol or 'UNK' for atom in lmp_types.values()]
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
        string += '%i bond types\n' % len(bond_types)
        string += '%i angle types\n' % len(angle_types)
        string += '%i dihedral types\n' % len(dihedral_types)
        string += '%i improper types\n' % len(improper_types)

        if not top.cell.is_rectangular:
            raise Exception('Triclinic box haven\'t been implemented')
        box = top.cell.get_size() * 10
        string += '\n0 %10.4f  xlo xhi\n' % box[0]
        string += '0 %10.4f  ylo yhi\n' % box[1]
        string += '0 %10.4f  zlo zhi\n' % box[2]

        string += '\nMasses\n\n'

        for i, (typ, atom) in enumerate(lmp_types.items()):
            string += '%4i %8.3f  # %8s\n' % (i + 1, atom.mass, typ)

        if len(bond_types) > 0:
            string += '\nBond Coeffs  # %s\n\n' % ('hybrid harmonic morse' if is_morse else 'harmonic')

        for i, bterm in enumerate(bond_types):
            bond_style = ''
            if is_morse:
                bond_style = 'morse' if type(bterm) is MorseBondTerm else 'harmonic'
            if type(bterm) is HarmonicBondTerm:
                string += '%4i %8s %12.6f %10.4f  # %s-%s\n' % (
                    i + 1, bond_style, bterm.k / 4.184 / 100, bterm.length * 10, bterm.type1, bterm.type2)
            elif type(bterm) is MorseBondTerm:
                alpha = (bterm.k / bterm.depth) ** 0.5
                string += '%4i %8s %12.6f %10.4f %10.4f  # %s-%s\n' % (
                    i + 1, bond_style, bterm.depth / 4.184, alpha / 10, bterm.length * 10, bterm.type1, bterm.type2)

        if len(angle_types) > 0:
            string += '\nAngle Coeffs  # %s\n\n' % ('hybrid harmonic sdk' if is_sdk else 'harmonic')

        for i, aterm in enumerate(angle_types):
            angle_style = ''
            if is_sdk:
                angle_style = 'sdk' if type(aterm) is SDKAngleTerm else 'harmonic'
            string += '%4i %8s %12.6f %10.4f  # %s-%s-%s\n' % (
                i + 1, angle_style, aterm.k / 4.184, aterm.theta * RAD2DEG, aterm.type1, aterm.type2, aterm.type3)

        if len(dihedral_types) > 0:
            string += '\nDihedral Coeffs  # opls\n\n'

        for i, dterm in enumerate(dihedral_types):
            if type(dterm) is PeriodicDihedralTerm:
                try:
                    dterm = dterm.to_opls_term()
                except:
                    raise Exception('Only OPLS type dihedral is implemented')
            k1, k2, k3, k4 = map(lambda x: x * 2 / 4.184, [dterm.k1, dterm.k2, dterm.k3, dterm.k4])
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
            typ = atom.type
            x, y, z = atom.position * 10
            string += '%8i %6i %6i %12.6f %10.4f %10.4f %10.4f  # %8s %8s\n' % (
                atom.id + 1, atom.molecule.id + 1, lmp_type_list.index(typ) + 1,
                atom.charge, x, y, z, atom.name, atom.molecule.name)

        if top.n_bond > 0:
            string += '\nBonds\n\n'

        for i, bond in enumerate(top.bonds):
            a1, a2 = bond.atom1, bond.atom2
            btype = bond_types.index(system.bond_terms[bond]) + 1
            string += '%6i %6i %6i %6i  # %s-%s\n' % (
                i + 1, btype, a1.id + 1, a2.id + 1, a1.name, a2.name)

        if top.n_angle > 0:
            string += '\nAngles\n\n'

        for i, angle in enumerate(top.angles):
            aterm = angle_types.index(system.angle_terms[angle]) + 1
            a1, a2, a3 = angle.atom1, angle.atom2, angle.atom3
            string += '%6i %6i %6i %6i %6i  # %s-%s-%s\n' % (
                i + 1, aterm, a1.id + 1, a2.id + 1, a3.id + 1, a1.name, a2.name, a3.name)

        if top.n_dihedral > 0:
            string += '\nDihedrals\n\n'

        for i, dihedral in enumerate(top.dihedrals):
            dtype = dihedral_types.index(system.dihedral_terms[dihedral]) + 1
            a1, a2, a3, a4 = dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4
            string += '%6i %6i %6i %6i %6i %6i  # %s-%s-%s-%s\n' % (
                i + 1, dtype, a1.id + 1, a2.id + 1, a3.id + 1, a4.id + 1,
                a1.name, a2.name, a3.name, a4.name)

        if top.n_improper > 0:
            string += '\nImpropers\n\n'

        for i, improper in enumerate(top.impropers):
            itype = improper_types.index(system.improper_terms[improper]) + 1
            a1, a2, a3, a4 = improper.atom1, improper.atom2, improper.atom3, improper.atom4
            string += '%6i %6i %6i %6i %6i %6i  # %s-%s-%s-%s\n' % (
                i + 1, itype, a2.id + 1, a3.id + 1, a1.id + 1, a4.id + 1,
                a2.name, a3.name, a1.name, a4.name)

        with open(data_out, 'wb') as f:
            f.write(string.encode())

        ### LAMMPS input script #################################################
        cmd_bond = 'hybrid harmonic morse' if is_morse else 'harmonic'
        cmd_angle = 'hybrid harmonic sdk' if is_sdk else 'harmonic'
        cmd_pair_base = 'lj/sdk' if is_sdk else 'lj/cut'
        cmd_pair_suffix = '/coul/long' if system.charged else ''
        cmd_mix = 'arithmetic' if system.ff.lj_mixing_rule == ForceField.LJ_MIXING_LB else 'geometric'
        cmd_tail = 'yes' if ff.vdw_long_range == ForceField.VDW_LONGRANGE_CORRECT else 'no'
        cmd_shift = 'yes' if ff.vdw_long_range == ForceField.VDW_LONGRANGE_SHIFT else 'no'
        cmd_comment_kspace = '# ' if not system.charged else ''

        cmd_paras = ''
        for i, (typ_i, atom) in enumerate(lmp_types.items()):
            for j, (typ_j, atom) in enumerate(lmp_types.items()):
                if j < i:
                    continue
                atype1 = lmp_types_atype[typ_i]
                atype2 = lmp_types_atype[typ_j]
                try:
                    # LAMMPS CG-SDK does not allow lj mixing. have to write it explicitly
                    vdw = ff.get_vdw_term(atype1, atype2, mixing=is_sdk)
                except:
                    continue

                if is_sdk:
                    if type(vdw) is MieTerm:
                        if vdw.is_sdk:
                            lj_style = f'lj{int(vdw.repulsion)}_{int(vdw.attraction)}'
                        else:
                            raise Exception('MieTerm other than SDK type are not implemented')
                    else:
                        lj_style = 'lj12_6'
                    cmd_paras += 'pair_coeff %4i %4i %8s %9.5f %8.4f  # %8s %8s %s\n' % (
                        i + 1, j + 1, lj_style, vdw.epsilon / 4.184, vdw.sigma * 10,
                        lmp_type_list[i], lmp_type_list[j], ','.join(vdw.comments))
                else:
                    cmd_paras += 'pair_coeff %4i %4i %9.5f %8.4f  # %8s %8s %s\n' % (
                        i + 1, j + 1, vdw.epsilon / 4.184, vdw.sigma * 10,
                        lmp_type_list[i], lmp_type_list[j], ','.join(vdw.comments))

        cmd_shake = 'fix SHAKE all shake 0.0001 20 0 b '
        if lmp_shake_bonds:
            cmd_shake += ' '.join(map(str, lmp_shake_bonds))
            if lmp_shake_angles:
                cmd_shake += ' a ' + ' '.join(map(str, lmp_shake_angles))
        else:
            cmd_shake = '# ' + cmd_shake
        timestep = 10.0 if is_sdk else 1.0

        string = f'''# created by mstk
units real
boundary p p p
atom_style full
bond_style {cmd_bond}
angle_style {cmd_angle}
dihedral_style opls
improper_style cvff
special_bonds lj 0 0 {ff.scale_14_vdw} coul 0 0 {ff.scale_14_coulomb}

pair_style {cmd_pair_base}{cmd_pair_suffix} {ff.vdw_cutoff * 10}
pair_modify mix {cmd_mix} tail {cmd_tail} shift {cmd_shift}
{cmd_comment_kspace}kspace_style pppm 1.0e-4

read_data {os.path.basename(data_out)}

{cmd_paras}

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
timestep {timestep}

fix NPT all npt temp $T $T {timestep * 100} iso $P $P {timestep * 500} fixedpoint 0 0 0

thermo_style custom step temp press pe emol evdwl v_elec density
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
