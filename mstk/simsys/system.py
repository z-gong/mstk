import copy
import itertools
from ..topology import *
from ..forcefield import *
from ..chem.constant import *
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
      You may need to call :func:`mstk.forcefield.ForceField.assign_charge`
      and :func:`mstk.forcefield.ForceField.assign_mass`
      before constructing the system to update the charges and masses of atoms in the topology.
    * If this is a Drude polarizable model, the Drude particles must be available in the topology,
      and the polarizability, thole screening and charges on Drude particles must have been correctly assigned.
      You may need to call :func:`mstk.topology.Topology.generate_drude_particles` to generate Drude particles.
    * If this is a virtual site model, the virtual site particles must be available in the topology,

    The provided force field must be able to fully describe the energy of the topology.
    If any term required by topology is not found in the force field, an Exception will be raised.

    Positions and unit cell are usually included in the topology.
    If not so, they can be provided by arguments `positions` and `cell`.
    The explicitly provided positions and cell will override the ones in topology.
    If positions are not included in the topology and not provided explicitly, an Exception will be raised.
    Unit cell is optional. If unit cell is not provided, cutoff will not be used for non-bonded interactions.

    If some bonded (bond, angle, dihedral, improper) parameters required by the topology are not found in the force field,
    `mstk` can try to transfer parameters from more general terms in the force field.
    e.g. if bond term between 'c_4o2' and 'n_3' is not found, the bond term between 'c_4' and 'n_3' will be used.
    This mechanism is initially designed for TEAM force field.
    It is disabled by default, and can be enabled by setting argument `transfer_bonded_terms` to True.

    Parameters
    ----------
    topology : Topology
    ff : ForceField
    positions : array_like, optional
    cell : UnitCell, optional
    ignore_missing_improper: bool
        If set to True, a zero-energy OplsImproperTerm will be created for missing improper terms in the FF
    transfer_bonded_terms: bool
    supress_pbc_warning: bool
        Set to True if you intend to build a vacuum system and don't want to hear warning about it.

    Attributes
    ----------
    charged : bool
        True if any atom in this system carries charge
    use_pbc : bool
        True if unit cell information provided and volume is not zero
    drude_pairs : dict, [Atom, Atom]
    vsite_pairs : dict, [Atom, Atom]
    constrain_bonds : dict, [Bond, float]
    constrain_angles : dict, [Angle, float]
    bond_terms : dict, [Bond, subclass of BondTerm]
    angle_terms : dict, [Angle, subclass of AngleTerm]
    dihedral_terms : dict, [Dihedral, subclass of DihedralTerm]
    improper_terms : dict, [Improper, subclass of ImproperTerm]
    polarizable_terms :  dict, [Atom, subclass of PolarizableTerm]
    vdw_classes : set of subclass of VdwTerm
    bond_classes : set of subclass of BondTerm
    angle_classes : set of subclass of AngleTerm
    dihedral_classes : set of subclass of DihedralTerm
    improper_classes : set of subclass of ImproperTerm
    polarizable_classes : set of subclass of PolarizableTerm
    ff_classes : set of subclass of FFTerm
    vsite_types : set of subclass of VirtualSite
    '''

    def __init__(self, topology, ff, positions=None, cell=None,
                 ignore_missing_improper=False,
                 transfer_bonded_terms=False,
                 suppress_pbc_warning=False):
        self._topology = copy.deepcopy(topology)
        self._ff = ForceField()

        if positions is not None:
            self._topology.set_positions(positions)
        elif not topology.has_position:
            raise Exception('Positions should be provided with topology or positions')

        if cell is not None:
            self._topology.cell = UnitCell(cell.vectors)

        if self._topology.n_atom == 0:
            logger.error('No atoms in the topology')
            raise Exception('Invalid topology')

        self.use_pbc = bool(self._topology.cell.volume != 0)  # convert numpy.bool_ to bool
        if not self.use_pbc and not suppress_pbc_warning:
            logger.warning('Periodic cell not provided or volume equal to zero')

        if all(atom.mass <= 0 for atom in topology.atoms):
            logger.error('All atoms have non-positive mass. '
                         'You may need to call ff.assign_mass(topology)')
            raise Exception('Invalid topology')

        self.charged = True
        if all(atom.charge == 0 for atom in topology.atoms):
            self.charged = False
            logger.warning('All atoms have zero charge. '
                           'You may need to call ff.assign_charge(topology)')

        charge_total = sum(atom.charge for atom in topology.atoms)
        if abs(charge_total) > 1E-8:
            logger.warning('System is not charge neutral. Charge = %.6f' % charge_total)

        # prebuild topology information
        self.drude_pairs: {Atom: Atom} = dict(self._topology.get_drude_pairs())  # {parent: atom_drude}
        self.vsite_pairs: {Atom: Atom} = dict(self._topology.get_virtual_site_pairs())  # {parent: atom_vsite}
        # records the types of VirtualSite in topology. It is NOT the VirtualSiteTerm in force field
        self.vsite_types = set([type(atom.virtual_site) for atom in self.vsite_pairs.values()])

        # bind FF terms with topological elements
        self.bond_terms: {Bond: BondTerm} = {}  # Drude bonds are not included
        self.angle_terms: {Angle: AngleTerm} = {}
        self.dihedral_terms: {Dihedral: DihedralTerm} = {}
        self.improper_terms: {Improper: ImproperTerm} = {}
        self.polarizable_terms: {Atom: PolarizableTerm} = {}  # key is parent Atom
        self.constrain_bonds: {Bond: float} = {}  # value is distance
        self.constrain_angles: {Angle: float} = {}  # value is 1-3 distance

        self.vdw_classes = set()
        self.bond_classes = set()
        self.angle_classes = set()
        self.dihedral_classes = set()
        self.improper_classes = set()
        self.polarizable_classes = set()
        self.ff_classes = set()

        self._extract_terms(ff, ignore_missing_improper, transfer_bonded_terms)

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

    def _extract_terms(self, ff, ignore_missing_improper=False, transfer_bonded_terms=False):
        '''
        Extract only the required terms from a force field.

        If any term required by topology is not found, an Exception will be raised.

        Parameters
        ----------
        ff : ForceField
        ignore_missing_improper : bool
        transfer_bonded_terms : bool
        '''
        from mstk.forcefield.dff_utils import dff_fuzzy_match

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
                self.bond_terms[bond] = bterm
                if bterm.fixed:
                    self.constrain_bonds[bond] = bterm.length

        _aterm_not_found = set()
        _aterm_transferred = set()
        for angle in self._topology.angles:
            _found = False
            ats_list = ff.get_eqt_for_angle(angle)
            _term = AngleTerm(*ats_list[0], 0)
            for ats in ats_list:
                try:
                    aterm = ff.get_angle_term(*ats)
                except FFTermNotFoundError:
                    pass
                else:
                    _found = True
                    break
            else:
                if transfer_bonded_terms:
                    aterm, score = dff_fuzzy_match(_term, ff)
                    if aterm is not None:
                        _found = True
                        _aterm_transferred.add((_term.name, aterm.name, score))
            if not _found:
                _aterm_not_found.add(_term.name)
            else:
                self._ff.add_term(aterm, replace=True)
                self.angle_terms[angle] = aterm
                if not aterm.fixed:
                    continue
                _bond1 = Bond(angle.atom1, angle.atom2)
                _bond2 = Bond(angle.atom2, angle.atom3)
                bond1 = next(b for b in angle.atom2.bonds if b.equals(_bond1))
                bond2 = next(b for b in angle.atom2.bonds if b.equals(_bond2))
                # constrain the angle only if two bonds are also constrained
                if bond1 in self.constrain_bonds and bond2 in self.constrain_bonds:
                    d1, d2 = self.constrain_bonds[bond1], self.constrain_bonds[bond2]
                    self.constrain_angles[angle] = math.sqrt(
                        d1 * d1 + d2 * d2 - 2 * d1 * d2 * math.cos(aterm.theta * PI / 180))

        _dterm_not_found = set()
        _dterm_transferred = set()
        _n_dihedral_to_remove = 0
        _dterm_to_remove = set()
        for mol in self._topology.molecules:
            _dihedral_to_remove = []
            _dihedrals_to_keep = []
            for dihedral in mol.dihedrals:
                _found = False
                ats_list = ff.get_eqt_for_dihedral(dihedral)
                _term = DihedralTerm(*ats_list[0])
                for ats in ats_list:
                    try:
                        dterm = ff.get_dihedral_term(*ats)
                    except FFTermNotFoundError:
                        pass
                    else:
                        _found = True
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
                        _dihedrals_to_keep.append(dihedral)
                        self._ff.add_term(dterm, replace=True)
                        self.dihedral_terms[dihedral] = dterm
            # extremely slow to remove lots of dihedrals by calling mol.remove_connectivity()
            mol._dihedrals = _dihedrals_to_keep

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
                except FFTermNotFoundError:
                    pass
                else:
                    _found = True
                    break
            else:
                if transfer_bonded_terms:
                    iterm, score = dff_fuzzy_match(_term, ff)
                    if iterm is not None:
                        _found = True
                        _iterm_transferred.add((_term.name, iterm.name, score))
            if not _found:
                if ignore_missing_improper:
                    iterm = OplsImproperTerm(*ats_list[0], 0.0)
                    _found = True
                else:
                    _iterm_not_found.add(_term.name)
            if _found:
                self._ff.add_term(iterm, replace=True)
                self.improper_terms[improper] = iterm

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

    def export_lammps(self, data_out='data.lmp', in_out='in.lmp', **kwargs):
        '''
        Generate input files for Lammps.

        Parameters
        ----------
        data_out : str
        in_out : str
        '''
        from .lmpexporter import LammpsExporter
        from .lmpdrudeexporter import LammpsDrudeExporter

        if DrudeTerm in self.ff_classes:
            LammpsDrudeExporter.export(self, data_out, in_out, **kwargs)
        else:
            LammpsExporter.export(self, data_out, in_out, **kwargs)

    def export_gromacs(self, gro_out='conf.gro', top_out='topol.top', mdp_out='grompp.mdp', **kwargs):
        '''
        Generate input files for GROMACS

        Parameters
        ----------
        gro_out : str or None
        top_out : str or None
        mdp_out : str or None
        '''
        from .gmxexporter import GromacsExporter
        GromacsExporter.export(self, gro_out, top_out, mdp_out, **kwargs)

    def export_namd(self, pdb_out='conf.pdb', psf_out='topol.psf', prm_out='ff.prm', **kwargs):
        '''
        Generate input files for NAMD

        Parameters
        ----------
        pdb_out : str or None
        psf_out : str or None
        prm_out : str or None
        '''
        from .namdexporter import NamdExporter
        NamdExporter.export(self, pdb_out, psf_out, prm_out, **kwargs)

    def to_omm_system(self, disable_inter_mol=False, **kwargs):
        '''
        Export this system to OpenMM system

        Parameters
        ----------
        disable_inter_mol : bool

        Returns
        -------
        omm_system : openmm.openmm.System
        '''
        from .ommexporter import OpenMMExporter
        return OpenMMExporter.export(self, disable_inter_mol=disable_inter_mol, **kwargs)

    def decompose_energy(self, verbose=True):
        '''
        Calculate the contribution of each energy term in this system

        Returns
        -------
        energies : Dict[str, float]
            The potential energies of different contributions
        '''
        try:
            from openmm import openmm as mm, app
            from mstk.ommhelper.utils import energy_decomposition
        except ImportError:
            raise Exception('OpenMM is required for energy decomposition')

        omm_sys = self.to_omm_system()
        integrator = mm.VerletIntegrator(0.001)
        platform = mm.Platform.getPlatformByName('CPU')
        sim = app.Simulation(self.topology.to_omm_topology(), omm_sys, integrator, platform)
        sim.context.setPositions(self.topology.positions)
        return energy_decomposition(sim, logger=logger if verbose else None)

    def minimize_energy(self, tolerance=100, verbose=True):
        '''
        Perform energy minimization on this system

        The position of each atom in the topology will be updated
        '''
        try:
            from openmm import openmm as mm, app
            from mstk.ommhelper.utils import minimize
        except ImportError:
            raise Exception('OpenMM is required for energy minimization')

        omm_sys = self.to_omm_system()
        integrator = mm.VerletIntegrator(0.001)
        platform = mm.Platform.getPlatformByName('CPU')
        sim = app.Simulation(self.topology.to_omm_topology(), omm_sys, integrator, platform)
        sim.context.setPositions(self.topology.positions)
        minimize(sim, tolerance=tolerance, logger=logger if verbose else None)
        positions = sim.context.getState(getPositions=True).getPositions(asNumpy=True)._value
        self.topology.set_positions(positions)

    def get_TIP4P_linear_coeffs(self, atom):
        '''
        Get the three linear coefficients to calculate the position of TIP4P virtual site from O, H and H positions

        Parameters
        ----------
        atom : Atom
            The virtual site atom in TIP4P water topology

        Returns
        -------
        params : tuple of float
        '''
        vsite: TIP4PSite = atom.virtual_site
        if type(vsite) is not TIP4PSite:
            raise Exception(f'{atom} is not a TIP4PSite')

        O, H1, H2 = vsite.parents
        mol = atom.molecule
        _bond_OH = Bond(O, H1)
        try:
            bond = next(b for b in mol.bonds if b.equals(_bond_OH))
        except StopIteration:
            raise Exception(f'Corresponding O-H bond not found for virtual atom {atom}')
        bterm = self.bond_terms[bond]

        _angle_HOH = Angle(H1, O, H2)
        try:
            angle = next(a for a in mol.angles if a.equals(_angle_HOH))
        except StopIteration:
            raise Exception(f'Corresponding H-O-H angle not found for virtual atom {atom}')
        aterm = self.angle_terms[angle]

        offset = atom.virtual_site.parameters[0]
        c = offset / bterm.length / 2 / math.cos(aterm.theta / 2 * PI / 180)

        return 1 - 2 * c, c, c
