import copy
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
    '''
    System is a combination of topology, forcefield and positions.
    When initiated, topology will be deep copied. And a compressed set of ff will be generated.
    Any modification to topology and ff after that will have no effect on the system.
    '''

    def __init__(self, topology: Topology, ff: FFSet, positions=None, cell=None):
        self._topology = copy.deepcopy(topology)
        self._ff = FFSet()

        if positions is not None:
            self._topology.set_positions(positions)
        elif not topology.has_position:
            raise Exception('Positions should be provided with topology or positions')

        if cell is not None:
            self._topology.cell = UnitCell(cell.vectors)
        if self._topology.cell.volume == 0:
            raise Exception('Periodic cell should be provided with topology or cell')

        self.extract_params(ff)
        self.charged = any(atom.charge != 0 for atom in self._topology.atoms)

    @property
    def topology(self):
        return self._topology

    @property
    def ff(self):
        return self._ff

    def extract_params(self, params: FFSet):
        self._ff.restore_settings(params.get_settings())

        self._bond_terms: {int: BondTerm} = {}  # key is id(Bond)
        self._angle_terms: {int: AngleTerm} = {}  # key is id(Angle)
        self._dihedral_terms: {int: DihedralTerm} = {}  # key is id(Dihedral)
        self._improper_terms: {int: ImproperTerm} = {}  # key is id(Improper)
        self._polarizable_terms: {Atom: PolarizableTerm} = {}  # key is parent Atom
        self._drude_pairs: {Atom: Atom} = dict(self._topology.get_drude_pairs())  # {parent: drude}
        self._constrain_bonds: {int: float} = {}  # key is id(Bond), value is distance
        self._constrain_angles: {int: float} = {}  # key is id(Angle), value is 1-3 distance

        for atom in self._topology.atoms:
            atype = params.atom_types.get(atom.type)
            if atype is None:
                raise Exception(f'AtomType {atom.type} not found in FF')

            if atype.name not in self._ff.atom_types:
                self._ff.add_term(atype)
                try:
                    vdw = params.get_vdw_term(atype, atype)
                except Exception as e:
                    raise Exception(f'Error assigning vdW parameters for {str(atom)}: {str(e)}')
                self._ff.add_term(vdw, ignore_duplicated=True)

        for atype1, atype2 in itertools.combinations(self._ff.atom_types.values(), 2):
            vdw = params.pairwise_vdw_terms.get(VdwTerm(atype1.eqt_vdw, atype2.eqt_vdw).name)
            if vdw is not None:
                self._ff.add_term(vdw, ignore_duplicated=True)

        for bond in self._topology.bonds:
            if bond.is_drude:
                # the Drude-parent bonds will be handled by DrudeForce below
                continue
            try:
                bterm = params.get_bond_term(bond)
            except Exception as e:
                raise Exception(f'Error assigning parameters for {str(bond)}: {str(e)}')
            self._ff.add_term(bterm, ignore_duplicated=True)
            self._bond_terms[id(bond)] = bterm
            if bterm.fixed:
                self._constrain_bonds[id(bond)] = bterm.length

        for angle in self._topology.angles:
            try:
                aterm = params.get_angle_term(angle)
            except Exception as e:
                raise Exception(f'Error assigning parameters for {str(angle)}: {str(e)}')
            self._ff.add_term(aterm, ignore_duplicated=True)
            self._angle_terms[id(angle)] = aterm
            if aterm.fixed:
                bond1 = next(b for b in self._topology.bonds if b == Bond(angle.atom1, angle.atom2))
                bond2 = next(b for b in self._topology.bonds if b == Bond(angle.atom2, angle.atom3))
                # constrain the angle only if two bonds are also constrained
                if id(bond1) in self._constrain_bonds and id(bond2) in self._constrain_bonds:
                    d1, d2 = self._constrain_bonds[id(bond1)], self._constrain_bonds[id(bond2)]
                    self._constrain_angles[id(angle)] = math.sqrt(
                        d1 * d1 + d2 * d2 - 2 * d1 * d2 * math.cos(aterm.theta * PI / 180))

        for dihedral in self._topology.dihedrals:
            try:
                dterm = params.get_dihedral_term(dihedral)
            except Exception as e:
                raise Exception(f'Error assign parameters for {str(dihedral)}: {str(e)}')
            self._ff.add_term(dterm, ignore_duplicated=True)
            self._dihedral_terms[id(dihedral)] = dterm

        for improper in self._topology.impropers:
            try:
                iterm = params.get_improper_term(improper)
            except Exception as e:
                raise Exception(f'Error assign parameters for {str(improper)}: {str(e)}')
            self._ff.add_term(iterm, ignore_duplicated=True)
            self._improper_terms[id(improper)] = iterm

        for parent in self._drude_pairs.keys():
            pterm = PolarizableTerm(params.atom_types[parent.type].eqt_polar)
            try:
                pterm = params.polarizable_terms[pterm.name]
            except:
                raise Exception(f'{str(pterm)} for {str(parent)} not found')
            self._ff.add_term(pterm, ignore_duplicated=True)
            self._polarizable_terms[parent] = pterm

        self.vdw_classes = {term.__class__ for term in self._ff.vdw_terms.values()}
        self.vdw_classes.update({term.__class__ for term in self._ff.pairwise_vdw_terms.values()})
        self.bond_classes = {term.__class__ for term in self._bond_terms.values()}
        self.angle_classes = {term.__class__ for term in self._angle_terms.values()}
        self.dihedral_classes = {term.__class__ for term in self._dihedral_terms.values()}
        self.improper_classes = {term.__class__ for term in self._improper_terms.values()}
        self.polarizable_classes = {term.__class__ for term in self._polarizable_terms.values()}

        self.ff_classes = (self.vdw_classes
                           .union(self.bond_classes)
                           .union(self.angle_classes)
                           .union(self.dihedral_classes)
                           .union(self.improper_classes)
                           .union(self.polarizable_classes))

    def export_lmp(self):
        from .lmpexporter import LammpsExporter
        LammpsExporter.export(self, data_out='data.lmp', in_out='in.lmp')

    def export_gmx(self):
        from .gmxexporter import GromacsExporter
        GromacsExporter.export(self, gro_out='conf.gro', top_out='topol.top', mdp_out='run.mdp')

    def to_omm_system(self):
        from .ommexporter import OpenMMExporter
        return OpenMMExporter.export(self)
