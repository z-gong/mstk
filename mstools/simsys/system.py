import copy
import itertools
import numpy as np
import sys
from ..topology import Topology, Atom, UnitCell, Psf, Bond, Angle, Dihedral, Improper
from ..trajectory import Frame, Trajectory, Gro
from ..forcefield import ForceField
from ..forcefield.ffterm import *
from ..constant import *
from .. import logger


class System():
    '''
    System is a combination of topology, forcefield and positions.
    When initiated, topology will be deep copied. And a compressed set of ff will be generated.
    Any modification to topology and ff after that will have no effect on the system.
    '''

    def __init__(self, topology: Topology, ff: ForceField, positions=None, cell=None):
        self._topology = copy.deepcopy(topology)
        self._ff = ForceField()

        if positions is not None:
            self._topology.set_positions(positions)
        elif not topology.has_position:
            raise Exception('Positions should be provided with topology or positions')

        if cell is not None:
            self._topology.cell = UnitCell(cell.vectors)
        if self._topology.cell.volume == 0:
            raise Exception('Periodic cell should be provided with topology or cell')

        if all(atom.mass <= 0 for atom in topology.atoms):
            logger.error('All atoms have non-positive mass. '
                         'You may need to call topology.assign_mass_from_ff(ff)')
            sys.exit(1)

        if all(atom.charge == 0 for atom in topology.atoms):
            self.charged = False
            logger.warning('All atoms have zero charge. '
                           'You may need to call topology.assign_charge_from_ff(ff)')
        else:
            self.charged = True

        self.extract_params(ff)

    @property
    def topology(self):
        return self._topology

    @property
    def ff(self):
        return self._ff

    def extract_params(self, params: ForceField):
        self._ff.restore_settings(params.get_settings())

        self._bond_terms: {int: BondTerm} = {}  # key is id(Bond), Drude bonds are not included
        self._angle_terms: {int: AngleTerm} = {}  # key is id(Angle)
        self._dihedral_terms: {int: DihedralTerm} = {}  # key is id(Dihedral)
        self._improper_terms: {int: ImproperTerm} = {}  # key is id(Improper)
        self._polarizable_terms: {Atom: PolarizableTerm} = {}  # key is parent Atom
        self._drude_pairs: {Atom: Atom} = dict(self._topology.get_drude_pairs())  # {parent: drude}
        self._constrain_bonds: {int: float} = {}  # key is id(Bond), value is distance
        self._constrain_angles: {int: float} = {}  # key is id(Angle), value is 1-3 distance

        _atype_not_found = set()
        for atom in self._topology.atoms:
            try:
                atype = params.atom_types[atom.type]
            except:
                _atype_not_found.add(atom.type)
            else:
                if atype.name not in self._ff.atom_types:
                    self._ff.add_term(atype)
        if _atype_not_found != set():
            logger.error('%i atom types not found in FF: %s' % (
                len(_atype_not_found), ' '.join(_atype_not_found)))
            raise Exception('Parameters missing')

        _vdw_not_found = set()
        for atype in self._ff.atom_types.values():
            try:
                vdw = params.get_vdw_term(atype, atype)
            except Exception as e:
                _vdw_not_found.add(atype.eqt_vdw)
            else:
                self._ff.add_term(vdw, replace=True)

        for atype1, atype2 in itertools.combinations(self._ff.atom_types.values(), 2):
            try:
                vdw = params.get_vdw_term(atype1, atype2, mixing=False)
            except:
                pass
            else:
                self._ff.add_term(vdw, replace=True)

        _bterm_not_found = set()
        for bond in self._topology.bonds:
            if bond.is_drude:
                # the Drude-parent bonds will be handled by DrudeForce below
                continue
            ats = params.get_eqt_for_bond(bond)
            try:
                bterm = params.get_bond_term(*ats)
            except:
                _bterm_not_found.add(BondTerm(*ats, 0).name)
            else:
                self._ff.add_term(bterm, replace=True)
                self._bond_terms[id(bond)] = bterm
                if bterm.fixed:
                    self._constrain_bonds[id(bond)] = bterm.length

        _aterm_not_found = set()
        for angle in self._topology.angles:
            ats = params.get_eqt_for_angle(angle)
            try:
                aterm = params.get_angle_term(*ats)
            except Exception as e:
                _aterm_not_found.add(AngleTerm(*ats, 0).name)
            else:
                self._ff.add_term(aterm, replace=True)
                self._angle_terms[id(angle)] = aterm
                if not aterm.fixed:
                    continue
                bond1 = next(b for b in angle.atom2.bonds if b == Bond(angle.atom1, angle.atom2))
                bond2 = next(b for b in angle.atom2.bonds if b == Bond(angle.atom2, angle.atom3))
                # constrain the angle only if two bonds are also constrained
                if id(bond1) in self._constrain_bonds and id(bond2) in self._constrain_bonds:
                    d1, d2 = self._constrain_bonds[id(bond1)], self._constrain_bonds[id(bond2)]
                    self._constrain_angles[id(angle)] = math.sqrt(
                        d1 * d1 + d2 * d2 - 2 * d1 * d2 * math.cos(aterm.theta * PI / 180))

        _dterm_not_found = set()
        for dihedral in self._topology.dihedrals:
            ats_list = params.get_eqt_for_dihedral(dihedral)
            for ats in ats_list:
                try:
                    dterm = params.get_dihedral_term(*ats)
                except Exception as e:
                    pass
                else:
                    self._ff.add_term(dterm, replace=True)
                    self._dihedral_terms[id(dihedral)] = dterm
                    break
            else:
                _dterm_not_found.add(DihedralTerm(*ats_list[0]).name)

        _iterm_not_found = set()
        for improper in self._topology.impropers:
            ats_list = params.get_eqt_for_improper(improper)
            for ats in ats_list:
                try:
                    iterm = params.get_improper_term(*ats)
                except Exception as e:
                    pass
                else:
                    self._ff.add_term(iterm, replace=True)
                    self._improper_terms[id(improper)] = iterm
                    break
            else:
                _iterm_not_found.add(ImproperTerm(*ats_list[0]).name)

        _pterm_not_found = set()
        for parent in self._drude_pairs.keys():
            pterm = PolarizableTerm(params.atom_types[parent.type].eqt_polar)
            try:
                pterm = params.polarizable_terms[pterm.name]
            except:
                _pterm_not_found.add(pterm.name)
            else:
                self._ff.add_term(pterm, replace=True)
                self._polarizable_terms[parent] = pterm

        ### print all the missing terms
        _not_found = False
        if _vdw_not_found != set():
            _not_found = True
            logger.error('%i vdW parameters not found in FF\n'
                         '        %s' % (len(_vdw_not_found),
                                         '\n        '.join(_vdw_not_found)))
        if _bterm_not_found != set():
            _not_found = True
            logger.error('%i bond parameters not found in FF\n'
                         '        %s' % (len(_bterm_not_found),
                                         '\n        '.join(_bterm_not_found)))
        if _aterm_not_found != set():
            _not_found = True
            logger.error('%i angle parameters not found in FF\n'
                         '        %s' % (len(_aterm_not_found),
                                         '\n        '.join(_aterm_not_found)))
        if _dterm_not_found != set():
            _not_found = True
            logger.error('%i dihedral parameters not found in FF\n'
                         '        %s' % (len(_dterm_not_found),
                                         '\n        '.join(_dterm_not_found)))
        if _iterm_not_found != set():
            _not_found = True
            logger.error('%i improper parameters not found in FF\n'
                         '        %s' % (len(_iterm_not_found),
                                         '\n        '.join(_iterm_not_found)))
        if _pterm_not_found != set():
            _not_found = True
            logger.error('%i polarizable parameters not found in FF\n'
                         '        %s' % (len(_pterm_not_found),
                                         '\n        '.join(_pterm_not_found)))

        if _not_found:
            raise Exception('Parameter missing')

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

    def export_lammps(self, data_out='data.lmp', in_out='in.lmp'):
        from .lmpexporter import LammpsExporter
        LammpsExporter.export(self, data_out, in_out)

    def export_gromacs(self, gro_out='conf.gro', top_out='topol.top', mdp_out='grompp.mdp'):
        from .gmxexporter import GromacsExporter
        GromacsExporter.export(self, gro_out, top_out, mdp_out)

    def export_charmm(self, pdb_out='conf.pdb', psf_out='topol.psf', prm_out='ff.prm'):
        from .charmmexporter import CharmmExporter
        CharmmExporter.export(self, pdb_out, psf_out, prm_out)

    def to_omm_system(self):
        from .ommexporter import OpenMMExporter
        return OpenMMExporter.export(self)
