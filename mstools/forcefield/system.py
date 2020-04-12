from ..topology import Topology, Atom
from ..trajectory import Frame, Trajectory
from .parameterset import ParameterSet
from .ffterm import *

try:
    from simtk import openmm as mm
    from simtk.openmm import app
except ImportError:
    OMM_EXIST = False
else:
    OMM_EXIST = True


class CONST:
    PI = 3.1415926535


class System():
    def __init__(self, topology: Topology, params: ParameterSet, frame: Frame):
        self._topology = topology
        self._params = params
        self._frame = frame

    def write_lmp(self):
        pass

    def write_gmx(self):
        pass

    def write_omm(self):
        pass

    def write_gro(self, file):
        from ..trajectory import Gro
        gro = Gro(file, 'w')
        gro.write_frame(self._topology, self._frame)
        gro.close()

    def write_psf(self, file):
        from ..topology import Psf
        psf = Psf(file, 'w')
        psf.is_drude = self._topology.is_drude
        psf.init_from_molecules(self._topology.molecules)
        psf.write()

    def to_omm_system(self):
        if not OMM_EXIST:
            raise Exception('cannot import openmm')

        omm_system = mm.System()
        for atom in self._topology.atoms:
            omm_system.addParticle(atom.mass)

        nbforce = mm.NonbondedForce()
        try:
            nbforce.setExceptionsUsePeriodicBoundaryConditions(True)
        except:
            print('warning: cannot apply periodic boundary conditions for exceptions. '
                  'be careful if there\'s infinite structures in the system')
        omm_system.addForce(nbforce)
        for atom in self._topology.atoms:
            atype = self._params.atom_types[atom.type]
            # use CustomNonbondedForce to handle vdW interactions
            nbforce.addParticle(atype.charge, 1.0, 0.0)
            nbforce.setUseDispersionCorrection(False)

        for bterm in self._params.bond_terms:
            if type(bterm) != HarmonicBondTerm:
                raise Exception('Bond terms other that HarmonicBondTerm haven\'t been implemented')
        bforce = mm.HarmonicBondForce()
        bforce.setUsesPeriodicBoundaryConditions(True)
        omm_system.addForce(bforce)
        for bond in self._topology.bonds:
            atype1: AtomType = self._params.atom_types[bond.atom1.type]
            atype2: AtomType = self._params.atom_types[bond.atom2.type]
            bterm = BondTerm(atype1.eqt_bond.name, atype2.eqt_bond.name, 0.0)
            bterm = self._params.bond_terms[bterm.name]
            if bterm.fixed:
                bforce.addBond(bond.atom1.id, bond.atom2.id, bterm.length, bterm.k)
            else:
                omm_system.addConstraint(bond.atom1, bond.atom2, bterm.length)

        for aterm in self._params.angle_terms:
            if type(aterm) != HarmonicAngleTerm:
                raise Exception(
                    'Angle terms other that HarmonicAngleTerm haven\'t been implemented')
        aforce = mm.HarmonicAngleForce()
        aforce.setUsesPeriodicBoundaryConditions(True)
        omm_system.addForce(aforce)
        for angle in self._topology.angles:
            atype1: AtomType = self._params.atom_types[angle.atom1.type]
            atype2: AtomType = self._params.atom_types[angle.atom2.type]
            atype3: AtomType = self._params.atom_types[angle.atom3.type]
            aterm = AngleTerm(atype1.eqt_angle_side.name, atype2.eqt_angle_center.name,
                              atype3.eqt_angle_side.name, 0.0)
            aterm = self._params.angle_terms[aterm.name]
            aforce.addAngle(angle.atom1.id, angle.atom2.id, angle.atom3.id,
                            aterm.theta * CONST.PI / 180, aterm.k)

        for dterm in self._params.dihedral_terms:
            if type(dterm) not in (PeriodicDihedralTerm, FourierDihedralTerm):
                raise Exception(
                    'Dihedral terms other that PeriodicDihedralTerm and FourierDihedralTerm'
                    'haven\'t been implemented')
        dforce = mm.PeriodicTorsionForce()
        dforce.setUsesPeriodicBoundaryConditions(True)
        omm_system.addForce(dforce)
        for dihedral in self._topology.dihedrals:
            dtype1: AtomType = self._params.atom_types[dihedral.atom1.type]
            dtype2: AtomType = self._params.atom_types[dihedral.atom2.type]
            dtype3: AtomType = self._params.atom_types[dihedral.atom3.type]
            dtype4: AtomType = self._params.atom_types[dihedral.atom4.type]
            dterm = DihedralTerm(dtype1.eqt_dihedral_side.name, dtype2.eqt_dihedral_center.name,
                                 dtype3.eqt_dihedral_center.name, dtype4.eqt_dihedral_side.name)

            dterm = self._params.dihedral_terms[dterm.name]
            ia1, ia2, ia3, ia4 = dihedral.atom1.id, dihedral.atom2.id, dihedral.atom3.id, dihedral.atom4.id
            if type(dterm) == PeriodicDihedralTerm:
                for para in dterm.parameters:
                    dforce.addTorsion(ia1, ia2, ia3, ia4, para.multiplicity,
                                      para.phi * CONST.PI / 180, para.k)
            elif type(dterm) == FourierDihedralTerm:
                if dterm.k1 != 0:
                    dforce.addTorsion(ia1, ia2, ia3, ia4, 1, 0, dterm.k1)
                if dterm.k2 != 0:
                    dforce.addTorsion(ia1, ia2, ia3, ia4, 2, CONST.PI, dterm.k2)
                if dterm.k3 != 0:
                    dforce.addTorsion(ia1, ia2, ia3, ia4, 3, 0, dterm.k3)
                if dterm.k4 != 0:
                    dforce.addTorsion(ia1, ia2, ia3, ia4, 4, CONST.PI, dterm.k4)

        for iterm in self._params.improper_terms:
            if type(iterm) not in (PeriodicImproperTerm,):
                raise Exception(
                    'Improper terms other that PeriodicImproperTerm haven\'t been implemented')
        iforce = mm.CustomTorsionForce('0.5*k*(1-cos(2*theta))')
        iforce.addPerTorsionParameter('k')
        iforce.setUsesPeriodicBoundaryConditions(True)
        omm_system.addForce(iforce)
        for improper in self._topology.impropers:
            itype1: AtomType = self._params.atom_types[improper.atom1.type]
            itype2: AtomType = self._params.atom_types[improper.atom2.type]
            itype3: AtomType = self._params.atom_types[improper.atom3.type]
            itype4: AtomType = self._params.atom_types[improper.atom4.type]
            iterm = ImproperTerm(itype1.eqt_improper_center.name, itype2.eqt_improper_side.name,
                                 itype3.eqt_improper_side.name, itype4.eqt_improper_side.name)
            iterm = self._params.improper_terms[iterm.name]
            iforce.addTorsion(improper.atom1.id, improper.atom2.id, improper.atom3.id,
                              improper.atom4.id, [iterm.k])

        cforce = mm.CustomNonbondedForce()
        cforce.setUseLongRangeCorrection(True)
        omm_system.addForce(cforce)

        pair12, pair13, pair14 = self._topology.get_12_13_14_exclusions()
        for atom1, atom2 in pair12 + pair13:
            cforce.addExclusion(atom1.id, atom2.id)
            nbforce.addException(atom1.id, atom2.id, 0.0, 1.0, 0.0)
        for atom1, atom2 in pair14:
            if self._params.scale_14_vdw == 1 and self._params.scale_14_coulomb == 1:
                continue
            if self._params.scale_14_vdw == 0 and self._params.scale_14_coulomb == 0:
                cforce.addExclusion(atom1.id, atom2.id)
                nbforce.addException(atom1.id, atom2.id, 0.0, 1.0, 0.0)
            else:
                charge_product = atom1.charge * atom2.charge * self._params.scale_14_coulomb
                if self._params.scale_14_vdw == 1:
                    nbforce.addException(atom1.id, atom2.id, charge_product, 1.0, 0.0)
                else:
                    cforce.addExclusion(atom1.id, atom2.id)
                    vdw = self._params.get_vdw_term(self._params.atom_types[atom1.type],
                                                    self._params.atom_types[atom2.type])
                    if type(vdw) == LJ126Term:
                        nbforce.addException(atom1.id, atom2.id, charge_product, vdw.sigma,
                                             vdw.epsilon * self._params.scale_14_vdw)
                    else:
                        raise Exception('vdW 1-4 scaling for non-LJ126 haven\'t been implemented')

        return omm_system
