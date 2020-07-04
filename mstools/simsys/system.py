import copy
import itertools
import sys
from ..topology import *
from ..forcefield import *
from ..constant import *
from .. import logger


class System():
    '''
    System is a combination of topology and forcefield, ready for simulation.

    When a system is being initiated, the topology is deep copied.
    A compressed set of force field is generated, which only contains the terms required by the topology.
    Any modification to the topology and force field after that will have no effect on the system.

    During the initialization, all information in the topology will be kept untouched, which means:

    * The atom type of all atoms in the topology must have been correctly assigned,
      so that the force field terms can be matched.
    * The charges and masses of all atoms in the topology must have been assigned.
      You may need to call :func:`mstools.topology.Topology.assign_charge_from_ff`
      and :func:`mstools.topology.Topology.assign_mass_from_ff`
      before constructing the system to update the charges and masses of atoms in the topology.
    * If this is a Drude polarizable model, the Drude particles must available in the topology,
      and the polarizability, thole screening and charges on Drude particles must have been correctly assigned.
      You may need to call :func:`mstools.topology.Topology.generate_drude_particles` to generate Drude particles.

    The provided force field must be able to fully describe the energy of the topology.
    If any term required by topology is not found in the force field, an Exception will be raised.

    Positions and unit cell are also required for construct a system.
    Usually they are included in the topology.
    If not so, they can be provided by arguments `positions` and `cell`.
    The explicitly provided positions and cell will override the ones in topology.
    If positions or unit cell is not included in the topology and not provided explicitly,
    an Exception will be raised.

    Parameters
    ----------
    topology : Topology
    ff : ForceField
    positions : array_like, optional
    cell : UnitCell, optional

    Attributes
    ----------
    charged : bool
        True if any atom in this system carries charge
    drude_pairs : dict, [Atom, Atom]
    constrain_bonds : dict, [int, float]
    constrain_angles : dict, [int, float]
    bond_terms : dict, [int, subclass of BondTerm]
    angle_terms : dict, [int, subclass of AngleTerm]
    dihedral_terms : dict, [int, subclass of DihedralTerm]
    improper_terms : dict, [int, subclass of ImproperTerm]
    polarizable_terms :  dict, [int, subclass of PolarizableTerm]
    vdw_classes : set of subclass of VdwTerm
    bond_classes : set of subclass of BondTerm
    angle_classes : set of subclass of AngleTerm
    dihedral_classes : set of subclass of DihedralTerm
    improper_classes : set of subclass of ImproperTerm
    polarizable_classes : set of subclass of PolarizableTerm
    ff_classes : set of subclass of FFTerm
    '''

    def __init__(self, topology, ff, positions=None, cell=None):
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

        self.charged = True
        if all(atom.charge == 0 for atom in topology.atoms):
            self.charged = False
            logger.warning('All atoms have zero charge. '
                           'You may need to call topology.assign_charge_from_ff(ff)')

        charge_total = sum(atom.charge for atom in topology.atoms)
        if abs(charge_total) > 1E-8:
            logger.warning('System is not charge neutral. Charge = %.6f' % charge_total)

        # bind FF terms with topological elements
        self.bond_terms: {int: BondTerm} = {}  # key is id(Bond), Drude bonds are not included
        self.angle_terms: {int: AngleTerm} = {}  # key is id(Angle)
        self.dihedral_terms: {int: DihedralTerm} = {}  # key is id(Dihedral)
        self.improper_terms: {int: ImproperTerm} = {}  # key is id(Improper)
        self.polarizable_terms: {Atom: PolarizableTerm} = {}  # key is parent Atom
        self.drude_pairs: {Atom: Atom} = dict(self._topology.get_drude_pairs())  # {parent: drude}
        self.constrain_bonds: {int: float} = {}  # key is id(Bond), value is distance
        self.constrain_angles: {int: float} = {}  # key is id(Angle), value is 1-3 distance

        self.vdw_classes = set()
        self.bond_classes = set()
        self.angle_classes = set()
        self.dihedral_classes = set()
        self.improper_classes = set()
        self.polarizable_classes = set()
        self.ff_classes = set()

        self._extract_terms(ff)

    @property
    def topology(self):
        '''
        The topology associated with the system

        Returns
        -------
        topology : Topology
        '''
        return self._topology

    @property
    def ff(self):
        '''
        The compressed force field associated with the system

        Returns
        -------
        ff : ForceField
        '''
        return self._ff

    def _extract_terms(self, ff: ForceField):
        '''
        Extract only the required terms from a force field.

        If any term required by topology is not found, an Exception will be raised.

        Parameters
        ----------
        ff : ForceField
        '''
        self._ff.restore_settings(ff.get_settings())

        _atype_not_found = set()
        for atom in self._topology.atoms:
            try:
                atype = ff.atom_types[atom.type]
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
                vdw = ff.get_vdw_term(atype, atype)
            except Exception as e:
                _vdw_not_found.add(atype.eqt_vdw)
            else:
                self._ff.add_term(vdw, replace=True)

        for atype1, atype2 in itertools.combinations(self._ff.atom_types.values(), 2):
            try:
                vdw = ff.get_vdw_term(atype1, atype2, mixing=False)
            except:
                pass
            else:
                self._ff.add_term(vdw, replace=True)

        _bterm_not_found = set()
        for bond in self._topology.bonds:
            if bond.is_drude:
                # the Drude-parent bonds will be handled by DrudeForce below
                continue
            ats = ff.get_eqt_for_bond(bond)
            try:
                bterm = ff.get_bond_term(*ats)
            except:
                _bterm_not_found.add(BondTerm(*ats, 0).name)
            else:
                self._ff.add_term(bterm, replace=True)
                self.bond_terms[id(bond)] = bterm
                if bterm.fixed:
                    self.constrain_bonds[id(bond)] = bterm.length

        _aterm_not_found = set()
        for angle in self._topology.angles:
            ats = ff.get_eqt_for_angle(angle)
            try:
                aterm = ff.get_angle_term(*ats)
            except Exception as e:
                _aterm_not_found.add(AngleTerm(*ats, 0).name)
            else:
                self._ff.add_term(aterm, replace=True)
                self.angle_terms[id(angle)] = aterm
                if not aterm.fixed:
                    continue
                bond1 = next(b for b in angle.atom2.bonds if b == Bond(angle.atom1, angle.atom2))
                bond2 = next(b for b in angle.atom2.bonds if b == Bond(angle.atom2, angle.atom3))
                # constrain the angle only if two bonds are also constrained
                if id(bond1) in self.constrain_bonds and id(bond2) in self.constrain_bonds:
                    d1, d2 = self.constrain_bonds[id(bond1)], self.constrain_bonds[id(bond2)]
                    self.constrain_angles[id(angle)] = math.sqrt(
                        d1 * d1 + d2 * d2 - 2 * d1 * d2 * math.cos(aterm.theta * PI / 180))

        _dterm_not_found = set()
        for dihedral in self._topology.dihedrals:
            ats_list = ff.get_eqt_for_dihedral(dihedral)
            for ats in ats_list:
                try:
                    dterm = ff.get_dihedral_term(*ats)
                except Exception as e:
                    pass
                else:
                    self._ff.add_term(dterm, replace=True)
                    self.dihedral_terms[id(dihedral)] = dterm
                    break
            else:
                _dterm_not_found.add(DihedralTerm(*ats_list[0]).name)

        _iterm_not_found = set()
        for improper in self._topology.impropers:
            ats_list = ff.get_eqt_for_improper(improper)
            for ats in ats_list:
                try:
                    iterm = ff.get_improper_term(*ats)
                except Exception as e:
                    pass
                else:
                    self._ff.add_term(iterm, replace=True)
                    self.improper_terms[id(improper)] = iterm
                    break
            else:
                _iterm_not_found.add(ImproperTerm(*ats_list[0]).name)

        _pterm_not_found = set()
        for parent in self.drude_pairs.keys():
            pterm = PolarizableTerm(ff.atom_types[parent.type].eqt_polar)
            try:
                pterm = ff.polarizable_terms[pterm.name]
            except:
                _pterm_not_found.add(pterm.name)
            else:
                self._ff.add_term(pterm, replace=True)
                self.polarizable_terms[parent] = pterm

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
        self.bond_classes = {term.__class__ for term in self.bond_terms.values()}
        self.angle_classes = {term.__class__ for term in self.angle_terms.values()}
        self.dihedral_classes = {term.__class__ for term in self.dihedral_terms.values()}
        self.improper_classes = {term.__class__ for term in self.improper_terms.values()}
        self.polarizable_classes = {term.__class__ for term in self.polarizable_terms.values()}
        self.ff_classes = (self.vdw_classes
                           .union(self.bond_classes)
                           .union(self.angle_classes)
                           .union(self.dihedral_classes)
                           .union(self.improper_classes)
                           .union(self.polarizable_classes))

    def export_lammps(self, data_out='data.lmp', in_out='in.lmp'):
        '''
        Generate input files for Lammps.

        Parameters
        ----------
        data_out : str
            Name of the data file
        in_out : str
            Name of the control script
        '''
        from .lmpexporter import LammpsExporter
        LammpsExporter.export(self, data_out, in_out)

    def export_gromacs(self, gro_out='conf.gro', top_out='topol.top', mdp_out='grompp.mdp'):
        '''
        Generate input files for GROMACS

        Parameters
        ----------
        gro_out : str or None
            Name of the GRO file
        top_out : str or None
            Name of the topol file
        mdp_out : str or None
            Name of the mdp control script
        '''
        from .gmxexporter import GromacsExporter
        GromacsExporter.export(self, gro_out, top_out, mdp_out)

    def export_charmm(self, pdb_out='conf.pdb', psf_out='topol.psf', prm_out='ff.prm'):
        '''
        Generate input files for CHARMM

        Parameters
        ----------
        pdb_out : str or None
            Name of the PDB file
        psf_out : str or None
            Name of the PSF file
        prm_out : str or None
            Name of the parameter file
        '''
        from .charmmexporter import CharmmExporter
        CharmmExporter.export(self, pdb_out, psf_out, prm_out)

    def to_omm_system(self):
        '''
        Export this system to OpenMM system

        Returns
        -------
        omm_system : simtk.openmm.System
        '''
        from .ommexporter import OpenMMExporter
        return OpenMMExporter.export(self)
