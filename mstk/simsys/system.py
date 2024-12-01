import itertools
from ..topology import *
from ..forcefield import *
from .. import logger


class System:
    '''
    System is a combination of topology and forcefield, ready for simulation.

    When a system being initialized, a compressed set of force field is constructed with the terms required by the topology.
    Note that the terms in the compressed force field are references to the terms in the original force field.
    A series of data structures are also constructed to map the topological elements to force field terms.
    Unless you are certain, DO NOT make any modification to the topology and force field after a system is initialized, as it may corrupt the system.

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
    * If this is a virtual site model, the virtual site particles must be available in the topology.

    The provided force field must be able to fully describe the energy of the topology.
    If any atom type from topology is not found in the FF, an Exception will be raised.
    If any vdw/bond/angle/dihedral/improper/polar term required by topology is not found in the FF,
    an Exception will be raised unless `allow_missing_terms` set to True.
    It won't check charge increment terms and virtual site terms, because they are already represented by the topology itself.

    Positions and unit cell should be included in the topology.
    If positions are not included in the topology, an Exception will be raised.
    Unit cell is optional. If unit cell is not provided, cutoff will not be used for non-bonded interactions.

    Parameters
    ----------
    topology : Topology
    ff : ForceField
    suppress_pbc_warning: bool
        Set to True if you intend to build a vacuum system and don't want to hear warning about it.
    allow_missing_terms: bool
        If set to True, won't raise Exception if FF terms other than AtomType are missing.
        This is useful if you want to check which FF terms are missing.

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
    atom_vdw_terms : dict, [Atom, subclass of VdwTerm]
    bond_terms : dict, [Bond, subclass of BondTerm]
    angle_terms : dict, [Angle, subclass of AngleTerm]
    dihedral_terms : dict, [Dihedral, subclass of DihedralTerm]
    improper_terms : dict, [Improper, subclass of ImproperTerm]
    polar_terms :  dict, [Atom, subclass of PolarTerm]
    missing_terms : list of FFTerm
    '''

    def __init__(self, topology, ff, suppress_pbc_warning=False, allow_missing_terms=False):
        self._topology = topology
        self._ff = ForceField()

        if not topology.has_position:
            raise Exception('Positions should be provided with topology')

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

        # map topological elements to FF terms
        self.atom_vdw_terms: {Atom: VdwTerm} = {}  # only same-type vdw terms are stored
        self.bond_terms: {Bond: BondTerm} = {}  # Drude bonds are not included
        self.angle_terms: {Angle: AngleTerm} = {}
        self.dihedral_terms: {Dihedral: DihedralTerm} = {}
        self.improper_terms: {Improper: ImproperTerm} = {}
        self.polar_terms: {Atom: PolarTerm} = {}  # key is parent Atom
        self.constrain_bonds: {Bond: float} = {}  # value is distance
        self.constrain_angles: {Angle: float} = {}  # value is 1-3 distance

        # missing FF terms
        self.missing_terms = []

        self._extract_terms(ff, allow_missing_terms)

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

    def _extract_terms(self, ff, allow_missing_terms=False):
        '''
        Extract only the required terms from a force field, and record the missing terms.

        Parameters
        ----------
        ff : ForceField
        allow_missing_terms : bool
        '''
        self._ff.restore_settings(ff.get_settings())

        _atype_not_found = set()
        for atom in self._topology.atoms:
            try:
                atype = ff.atom_types[atom.type]
            except:
                _atype_not_found.add(atom.type)
            else:
                self._ff.add_term(atype, replace=True)
        if _atype_not_found:
            logger.error('%i atom types not found in FF: %s' % (len(_atype_not_found), ' '.join(_atype_not_found)))
            raise Exception('Incomplete force field')

        _vdw_not_found = {}
        for atype in self._ff.atom_types.values():
            try:
                vdw = ff.get_vdw_term(atype, atype)
            except Exception as e:
                _vdw = VdwTerm(atype.eqt_vdw, atype.eqt_vdw)
                _vdw_not_found[_vdw.name] = _vdw
            else:
                self._ff.add_term(vdw, replace=True)

        for atype1, atype2 in itertools.combinations(self._ff.atom_types.values(), 2):
            try:
                vdw = ff.get_vdw_term(atype1, atype2, mixing=False)
            except:
                if self._ff.lj_mixing_rule == self._ff.LJ_MIXING_NONE:
                    # all pairwise vdW terms are explicitly required
                    _vdw = VdwTerm(atype1.eqt_vdw, atype2.eqt_vdw)
                    _vdw_not_found[_vdw.name] = _vdw
            else:
                self._ff.add_term(vdw, replace=True)

        # map atom to self vdW term
        for atom in self._topology.atoms:
            atype = self._ff.atom_types[atom.type]
            try:
                vdw = self._ff.get_vdw_term(atype, atype)
            except FFTermNotFoundError:
                pass
            else:
                self.atom_vdw_terms[atom] = vdw

        _bterm_not_found = {}
        for bond in self._topology.bonds:
            if bond.is_drude:
                # the Drude-parent bonds will be handled by DrudeForce below
                continue
            ats = ff.get_eqt_for_bond(bond)
            _term = BondTerm(*ats, 0)
            try:
                bterm = ff.get_bond_term(*ats)
            except FFTermNotFoundError:
                _bterm_not_found[_term.name] = _term
            else:
                self._ff.add_term(bterm, replace=True)
                self.bond_terms[bond] = bterm
                if bterm.fixed:
                    self.constrain_bonds[bond] = bterm.length

        _aterm_not_found = {}
        for angle in self._topology.angles:
            _found = False
            ats_list = ff.get_eqt_list_for_angle(angle)
            _term = AngleTerm(*ats_list[0], 0)
            for ats in ats_list:
                try:
                    aterm = ff.get_angle_term(*ats)
                except FFTermNotFoundError:
                    pass
                else:
                    _found = True
                    break
            if not _found:
                _aterm_not_found[_term.name] = _term
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
                    self.constrain_angles[angle] = math.sqrt(d1 * d1 + d2 * d2 - 2 * d1 * d2 * math.cos(aterm.theta))

        _dterm_not_found = {}
        for mol in self._topology.molecules:
            for dihedral in mol.dihedrals:
                _found = False
                ats_list = ff.get_eqt_list_for_dihedral(dihedral)
                _term = DihedralTerm(*ats_list[0])
                for ats in ats_list:
                    try:
                        dterm = ff.get_dihedral_term(*ats)
                    except FFTermNotFoundError:
                        pass
                    else:
                        _found = True
                        break
                if not _found:
                    _dterm_not_found[_term.name] = _term
                else:
                    self._ff.add_term(dterm, replace=True)
                    self.dihedral_terms[dihedral] = dterm

        _iterm_not_found = {}
        for improper in self._topology.impropers:
            _found = False
            ats_list = ff.get_eqt_list_for_improper(improper)
            _term = ImproperTerm(*ats_list[0])
            for ats in ats_list:
                try:
                    iterm = ff.get_improper_term(*ats)
                except FFTermNotFoundError:
                    pass
                else:
                    _found = True
                    break
            if not _found:
                _iterm_not_found[_term.name] = _term
            else:
                self._ff.add_term(iterm, replace=True)
                self.improper_terms[improper] = iterm

        _pterm_not_found = {}
        for parent in self.drude_pairs.keys():
            _term = PolarTerm(ff.atom_types[parent.type].eqt_polar)
            try:
                pterm = ff.polar_terms[_term.name]
            except:
                _pterm_not_found[_term.name] = _term
            else:
                self._ff.add_term(pterm, replace=True)
                self.polar_terms[parent] = pterm

        ### print all the missing terms
        for d in _vdw_not_found, _bterm_not_found, _aterm_not_found, _dterm_not_found, _iterm_not_found, _pterm_not_found:
            self.missing_terms.extend(d.values())

        if self.missing_terms:
            msg = f'{len(self.missing_terms)} terms not found in FF'
            for term in self.missing_terms:
                msg += f'\n        {term}'
            if allow_missing_terms:
                logger.warning(msg)
            else:
                logger.error(msg)
                raise Exception('Incomplete force field')

    def export_lammps(self, data_out='data.lmp', in_out='in.lmp', **kwargs):
        '''
        Generate input files for Lammps.

        Parameters
        ----------
        data_out : str
        in_out : str
        '''
        from .lmpexporter import LammpsExporter
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

    def export_namd(self, pdb_out='conf.pdb', psf_out='top.psf', prm_out='ff.prm', **kwargs):
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
            Disable inter-molecular interactions. This is useful for optimizing a batch of molecules

        Returns
        -------
        omm_system : openmm.openmm.System
        '''
        from .ommexporter import OpenMMExporter
        return OpenMMExporter.export(self, disable_inter_mol=disable_inter_mol, **kwargs)

    def decompose_energy(self, disable_inter_mol=False, verbose=True):
        '''
        Calculate the contribution of each energy term in this system

        Parameters
        ----------
        disable_inter_mol : bool
        verbose : bool

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

        omm_sys = self.to_omm_system(disable_inter_mol=disable_inter_mol)
        integrator = mm.VerletIntegrator(0.001)
        platform = mm.Platform.getPlatformByName('CPU')
        sim = app.Simulation(self.topology.to_omm_topology(), omm_sys, integrator, platform)
        sim.context.setPositions(self.topology.positions)
        return energy_decomposition(sim, logger=logger if verbose else None)

    def minimize_energy(self, disable_inter_mol=False, tolerance=100.0, verbose=True):
        '''
        Perform energy minimization on this system.

        The position of each atom in the topology will be updated.

        Parameters
        ----------
        disable_inter_mol : bool
        tolerance : float
        verbose : bool
        '''
        try:
            from openmm import openmm as mm, app
            from mstk.ommhelper.utils import minimize
        except ImportError:
            raise Exception('OpenMM is required for energy minimization')

        omm_sys = self.to_omm_system(disable_inter_mol=disable_inter_mol)
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
        c = offset / bterm.length / 2 / math.cos(aterm.theta / 2)

        return 1 - 2 * c, c, c
