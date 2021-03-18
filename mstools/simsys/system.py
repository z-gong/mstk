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

    If some bonded (bond, angle, dihedral, improper) parameters required by the topology are not found in the force field,
    `mstools` can try to transfer parameters from more general terms in the force field.
    e.g. if bond term between 'c_4o2' and 'n_3' is not found, the bond term between 'c_4' and 'n_3' will be used.
    This mechanism is initially designed for TEAM force field.
    It is disabled by default, and can be enabled by setting argument `transfer_bonded_terms` to True.

    Parameters
    ----------
    topology : Topology
    ff : ForceField
    positions : array_like, optional
    cell : UnitCell, optional
    transfer_bonded_terms: bool, optional
    supress_pbc_warning: bool, optional
        Set to True if you intend to build a vacuum system and don't want to hear warning about it.

    Attributes
    ----------
    charged : bool
        True if any atom in this system carries charge
    use_pbc : bool
        True if unit cell information provided and volume is not zero
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

    def __init__(self, topology, ff, positions=None, cell=None, transfer_bonded_terms=False,
                 suppress_pbc_warning=False):
        self._topology = copy.deepcopy(topology)
        self._ff = ForceField()

        if positions is not None:
            self._topology.set_positions(positions)
        elif not topology.has_position:
            raise Exception('Positions should be provided with topology or positions')

        if cell is not None:
            self._topology.cell = UnitCell(cell.vectors)

        self.use_pbc = bool(self._topology.cell.volume != 0)  # convert numpy.bool_ to bool
        if not self.use_pbc and not suppress_pbc_warning:
            logger.warning('Periodic cell not provided or volume equal to zero')

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

        # prebuild topology information
        self.drude_pairs: {Atom: Atom} = dict(self._topology.get_drude_pairs())  # {parent: atom_drude}
        self.vsite_pairs: {Atom: Atom} = dict(self._topology.get_virtual_site_pairs())  # {parent: atom_vsite}
        # records the types of VirtualSite in topology. It is NOT the VirtualSiteTerm in force field
        self.vsite_types = set([type(atom.virtual_site) for atom in self.vsite_pairs.values()])

        # bind FF terms with topological elements
        self.bond_terms: {int: BondTerm} = {}  # key is id(Bond), Drude bonds are not included
        self.angle_terms: {int: AngleTerm} = {}  # key is id(Angle)
        self.dihedral_terms: {int: DihedralTerm} = {}  # key is id(Dihedral)
        self.improper_terms: {int: ImproperTerm} = {}  # key is id(Improper)
        self.polarizable_terms: {Atom: PolarizableTerm} = {}  # key is parent Atom
        self.constrain_bonds: {int: float} = {}  # key is id(Bond), value is distance
        self.constrain_angles: {int: float} = {}  # key is id(Angle), value is 1-3 distance

        self.vdw_classes = set()
        self.bond_classes = set()
        self.angle_classes = set()
        self.dihedral_classes = set()
        self.improper_classes = set()
        self.polarizable_classes = set()
        self.ff_classes = set()

        self._extract_terms(ff, transfer_bonded_terms)

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

    def _extract_terms(self, ff, transfer_bonded_terms=False):
        '''
        Extract only the required terms from a force field.

        If any term required by topology is not found, an Exception will be raised.

        Parameters
        ----------
        ff : ForceField
        transfer_bonded_terms : bool
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
        _bterm_transferred = set()
        for bond in self._topology.bonds:
            if bond.is_drude:
                # the Drude-parent bonds will be handled by DrudeForce below
                continue
            _found = False
            ats = ff.get_eqt_for_bond(bond)
            _term = BondTerm(*ats, 0)
            try:
                bterm = ff.get_bond_term(*ats)
                _found = True
            except:
                if transfer_bonded_terms:
                    bterm, score = dff_fuzzy_match(_term, ff)
                    if bterm is not None:
                        _found = True
                        _bterm_transferred.add((_term.name, bterm.name, score))
            if not _found:
                _bterm_not_found.add(_term.name)
            else:
                self._ff.add_term(bterm, replace=True)
                self.bond_terms[id(bond)] = bterm
                if bterm.fixed:
                    self.constrain_bonds[id(bond)] = bterm.length

        _aterm_not_found = set()
        _aterm_transferred = set()
        for angle in self._topology.angles:
            _found = False
            ats = ff.get_eqt_for_angle(angle)
            _term = AngleTerm(*ats, 0)
            try:
                aterm = ff.get_angle_term(*ats)
                _found = True
            except:
                if transfer_bonded_terms:
                    aterm, score = dff_fuzzy_match(_term, ff)
                    if aterm is not None:
                        _found = True
                        _aterm_transferred.add((_term.name, aterm.name, score))
            if not _found:
                _aterm_not_found.add(_term.name)
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
        _dterm_transferred = set()
        _n_dihedral_to_remove = 0
        _dterm_to_remove = set()
        for mol in self._topology.molecules:
            _dihedral_to_remove = []
            for dihedral in mol.dihedrals:
                _found = False
                ats_list = ff.get_eqt_for_dihedral(dihedral)
                _term = DihedralTerm(*ats_list[0])
                for ats in ats_list:
                    try:
                        dterm = ff.get_dihedral_term(*ats)
                        _found = True
                    except:
                        pass
                    else:
                        break
                else:
                    if transfer_bonded_terms:
                        dterm, score = dff_fuzzy_match(_term, ff)
                        if dterm is not None:
                            _found = True
                            _dterm_transferred.add((_term.name, dterm.name, score))
                if not _found:
                    _dterm_not_found.add(_term.name)
                else:
                    if dterm.is_zero:
                        # remove dihedrals that have zero energy contributions
                        # in case linear groups like alkyne and nitrile give energy of NaN
                        _dihedral_to_remove.append(dihedral)
                        _n_dihedral_to_remove += 1
                        _dterm_to_remove.add(dterm.name)
                    else:
                        self._ff.add_term(dterm, replace=True)
                        self.dihedral_terms[id(dihedral)] = dterm
            # remove dihedrals from this molecule that have zero energy contributions
            for dihedral in _dihedral_to_remove:
                mol.remove_connectivity(dihedral)

        if _n_dihedral_to_remove > 0:
            logger.warning('%i dihedrals removed because they have zero energy contributions\n'
                           '        %s' % (_n_dihedral_to_remove,
                                           '\n        '.join(_dterm_to_remove)))

        _iterm_not_found = set()
        _iterm_transferred = set()
        for improper in self._topology.impropers:
            _found = False
            ats_list = ff.get_eqt_for_improper(improper)
            _term = ImproperTerm(*ats_list[0])
            for ats in ats_list:
                try:
                    iterm = ff.get_improper_term(*ats)
                    _found = True
                except:
                    pass
                else:
                    break
            else:
                if transfer_bonded_terms:
                    iterm, score = dff_fuzzy_match(_term, ff)
                    if iterm is not None:
                        _found = True
                        _iterm_transferred.add((_term.name, iterm.name, score))
            if not _found:
                _iterm_not_found.add(_term.name)
            else:
                self._ff.add_term(iterm, replace=True)
                self.improper_terms[id(improper)] = iterm

        _pterm_not_found = set()
        for parent in self.drude_pairs.keys():
            _term = PolarizableTerm(ff.atom_types[parent.type].eqt_polar)
            try:
                pterm = ff.polarizable_terms[_term.name]
            except:
                _pterm_not_found.add(_term.name)
            else:
                self._ff.add_term(pterm, replace=True)
                self.polarizable_terms[parent] = pterm

        ### print all transferred terms
        _n_transferred = 0
        msg = ''
        for k, v in {'Bond'    : _bterm_transferred,
                     'Angle'   : _aterm_transferred,
                     'Dihedral': _dterm_transferred,
                     'Improper': _iterm_transferred}.items():
            for i in v:
                _n_transferred += 1
                msg += '        %-10s %-20s %-20s %i\n' % (k, *i)
        if _n_transferred > 0:
            logger.warning(
                '%i bonded terms not found in FF. Transferred from similar terms with score\n' % _n_transferred + msg)

        ### print all the missing terms
        _n_missing = 0
        msg = ''
        for k, v in {'vdW'        : _vdw_not_found,
                     'Bond'       : _bterm_not_found,
                     'Angle'      : _aterm_not_found,
                     'Dihedral'   : _dterm_not_found,
                     'Improper'   : _iterm_not_found,
                     'Polarizable': _pterm_not_found}.items():
            for i in v:
                _n_missing += 1
                msg += '        %-10s %-20s\n' % (k, i)
        if _n_missing > 0:
            logger.error('%i terms not found in FF\n' % _n_missing + msg)
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
        in_out : str
        '''
        from .lmpexporter import LammpsExporter
        LammpsExporter.export(self, data_out, in_out)

    def export_gromacs(self, gro_out='conf.gro', top_out='topol.top', mdp_out='grompp.mdp'):
        '''
        Generate input files for GROMACS

        Parameters
        ----------
        gro_out : str or None
        top_out : str or None
        mdp_out : str or None
        '''
        from .gmxexporter import GromacsExporter
        GromacsExporter.export(self, gro_out, top_out, mdp_out)

    def export_charmm(self, pdb_out='conf.pdb', psf_out='topol.psf', prm_out='ff.prm'):
        '''
        Generate input files for CHARMM

        Parameters
        ----------
        pdb_out : str or None
        psf_out : str or None
        prm_out : str or None
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

    def get_TIP4P_linear_coeffs(self, atom: Atom):
        vsite: TIP4PSite = atom.virtual_site
        if type(vsite) != TIP4PSite:
            raise Exception(f'{repr(atom)} is not a TIP4PSite')

        O, H1, H2 = vsite.parents
        mol = atom.molecule
        bond_OH = Bond(O, H1)
        for bond in mol.bonds:
            if bond == bond_OH:
                bterm = self.bond_terms[id(bond)]
                break
        else:
            raise Exception(f'Corresponding O-H bond not found for virtual atom {repr(atom)}')

        angle_HOH = Angle(H1, O, H2)
        for angle in mol.angles:
            if angle == angle_HOH:
                aterm = self.angle_terms[id(angle)]
                break
        else:
            raise Exception(f'Corresponding H-O-H angle not found for virtual atom {repr(atom)}')

        offset = atom.virtual_site.parameters[0]
        c = offset / bterm.length / 2 / math.cos(aterm.theta / 2 * PI / 180)

        return 1 - 2 * c, c, c
