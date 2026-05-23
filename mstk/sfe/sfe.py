"""
Solvation free energy (SFE) via alchemical decoupling using mstk + OpenMM.

Decouples a single molecule from its environment in two phases:
  Phase 1: Coulomb (lambda_coul 1->0), LJ fully on
  Phase 2: LJ with soft-core (lambda_vdw 1->0), Coulomb already off
"""

import numpy as np
from openmm import openmm as mm

from ..forcefield import ForceField, LJ126Term
from ..chem import constant


class SFEManager:
    """
    Alchemical system for solvation free energy decoupling.

    Takes an mstk System and builds an OpenMM system with alchemical forces
    for decoupling a single molecule from its environment.
    Two-phase protocol: Coulomb first, then LJ with soft-core.

    Parameters
    ----------
    system : mstk.simsys.System
        The molecular system (topology + force field).
    mol_id : int
        Index of the molecule to decouple (0-based).
    n_lambda : int, optional
        Number of lambda windows (default 16).

    Attributes
    ----------
    omm_system : openmm.System
        The alchemical OpenMM system, ready to create a Simulation.
    schedule : list of (float, float)
        Lambda schedule. Each element is (lambda_coul, lambda_vdw).
    mol_atoms : list of int
        Global atom indices of the decoupled molecule.
    topology : mstk.topology.Topology
        Reference to the original topology.

    Examples
    --------
    >>> # Create alchemical manager
    >>> sfe = SFEManager(system, mol_id=0, n_lambda=16)
    >>>
    >>> # Create simulation with the alchemical OpenMM system
    >>> integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.ps, 0.002*unit.ps)
    >>> sim = app.Simulation(top.to_omm_topology(), sfe.omm_system, integrator)
    >>> sim.context.setPositions(top.positions)
    >>>
    >>> # Set lambda state for a specific window
    >>> sfe.set_lambda_state(sim.context, 5)
    >>>
    >>> # Evaluate dU to all windows at current configuration
    >>> U = []
    >>> for k in range(16):
    ...     sfe.set_lambda_state(sim.context, k)
    ...     U.append(sim.context.getState(getEnergy=True).getPotentialEnergy())
    >>>
    >>> # After collecting dU from all windows, run MBAR
    >>> mu_ex, error, dG_adj, dG_adj_err = SFEManager.compute_mbar(u_kn, N_k, 300.0)
    """

    def __init__(self, system, mol_id, n_lambda=16):
        top = system.topology
        ff = system.ff

        unsupported_vdw = ff.vdw_term_classes - {LJ126Term}
        if unsupported_vdw:
            raise ValueError(
                f"Unsupported VdW term classes: "
                f"{', '.join(c.__name__ for c in unsupported_vdw)}. "
                f"Alchemical soft-core only supports LJ126Term."
            )
        if ff.vdw_long_range == ForceField.VDW_LONGRANGE_SHIFT:
            raise ValueError(
                "vdw_long_range='shift' is not supported. "
                "Alchemical force uses unshifted LJ potential."
            )
        if ff.is_polarizable:
            raise ValueError(
                "Polarizable force fields (Drude) are not supported. "
                "Alchemical decoupling does not handle induced dipoles."
            )
        if ff.has_virtual_site:
            raise ValueError(
                "Virtual site force fields (e.g. TIP4P) are not supported. "
                "Alchemical decoupling does not handle virtual site charges."
            )

        self.topology = top
        self.schedule = self.default_lambda_schedule(n_lambda)
        mol = top.molecules[mol_id]

        self.mol_atoms = [atom.id for atom in mol.atoms]
        mol_atom_set = set(self.mol_atoms)
        all_atom_ids = set(range(top.n_atom))
        env_atom_ids = all_atom_ids - mol_atom_set

        self.omm_system = system.to_omm_system()

        # Find forces.
        self._nbforce = None
        cnb_forces = []
        for force in self.omm_system.getForces():
            if isinstance(force, mm.NonbondedForce):
                self._nbforce = force
            elif isinstance(force, mm.CustomNonbondedForce):
                cnb_forces.append(force)

        if self._nbforce is None:
            raise ValueError("No NonbondedForce found in the system")
        if not cnb_forces:
            raise ValueError("No CustomNonbondedForce found in the system")

        # Store original charges for M atoms.
        self._original_charges = {}
        for idx in self.mol_atoms:
            q, _, _ = self._nbforce.getParticleParameters(idx)
            self._original_charges[idx] = q._value

        # Identify beyond-1-4 pairs within M.
        distance_matrix = mol.get_distance_matrix(max_bond=3)
        beyond_14_pairs = []
        for i in range(mol.n_atom):
            for j in range(i + 1, mol.n_atom):
                if distance_matrix[i, j] == 0:
                    qi = mol.atoms[i].charge
                    qj = mol.atoms[j].charge
                    beyond_14_pairs.append((self.mol_atoms[i], self.mol_atoms[j], qi * qj))

        # --- Coulomb correction CustomBondForce ---
        coul_correction = mm.CustomBondForce(
            f'(1 - lambda_coul^2) * {constant.ONE_4PI_EPS0} * qq / r'
        )
        coul_correction.addGlobalParameter('lambda_coul', 1.0)
        coul_correction.addPerBondParameter('qq')
        coul_correction.setUsesPeriodicBoundaryConditions(True)
        coul_correction.setForceGroup(10)
        coul_correction.setName('CoulCorrection')
        for (ai, aj, qq) in beyond_14_pairs:
            coul_correction.addBond(ai, aj, [qq])
        self.omm_system.addForce(coul_correction)

        # --- LJ decoupling ---
        original_cnb = cnb_forces[0]
        cutoff = original_cnb.getCutoffDistance()
        n_type = None
        A_values = None
        B_values = None

        for func_idx in range(original_cnb.getNumTabulatedFunctions()):
            name = original_cnb.getTabulatedFunctionName(func_idx)
            func = original_cnb.getTabulatedFunction(func_idx)
            if isinstance(func, mm.Discrete2DFunction) and name == 'A':
                xsize, _, A_values = func.getFunctionParameters()
                n_type = xsize
            elif isinstance(func, mm.Discrete2DFunction) and name == 'B':
                _, _, B_values = func.getFunctionParameters()

        if n_type is None or A_values is None or B_values is None:
            raise ValueError("Cannot extract A/B tables from CustomNonbondedForce")

        # Compute S table (sigma^6 for soft-core denominator).
        S_values = [0.0] * (n_type * n_type)
        for idx in range(n_type * n_type):
            if B_values[idx] > 0:
                S_values[idx] = A_values[idx] / B_values[idx]

        # Restrict original force to (env,env) + (M,M).
        original_cnb.addInteractionGroup(env_atom_ids, env_atom_ids)
        original_cnb.addInteractionGroup(mol_atom_set, mol_atom_set)

        # Build alchemical LJ force (M-env, soft-core).
        alch_lj = mm.CustomNonbondedForce(
            'lambda_vdw * (A(type1,type2) * invR6_sc * invR6_sc'
            ' - B(type1,type2) * invR6_sc);'
            'invR6_sc = 1 / (alpha * (1-lambda_vdw)^2 * S(type1,type2) + r^6);'
            'alpha = 0.5'
        )
        alch_lj.addGlobalParameter('lambda_vdw', 1.0)
        alch_lj.addPerParticleParameter('type')
        alch_lj.addTabulatedFunction('A', mm.Discrete2DFunction(n_type, n_type, A_values))
        alch_lj.addTabulatedFunction('B', mm.Discrete2DFunction(n_type, n_type, B_values))
        alch_lj.addTabulatedFunction('S', mm.Discrete2DFunction(n_type, n_type, S_values))

        for i in range(original_cnb.getNumParticles()):
            alch_lj.addParticle(original_cnb.getParticleParameters(i))
        for i in range(original_cnb.getNumExclusions()):
            alch_lj.addExclusion(*original_cnb.getExclusionParticles(i))

        alch_lj.addInteractionGroup(mol_atom_set, env_atom_ids)
        alch_lj.setNonbondedMethod(original_cnb.getNonbondedMethod())
        alch_lj.setCutoffDistance(cutoff)
        alch_lj.setUseLongRangeCorrection(True)
        alch_lj.setForceGroup(11)
        alch_lj.setName('AlchLJ')
        self.omm_system.addForce(alch_lj)

    def set_lambda_state(self, context, window):
        """
        Set the alchemical state of the context by window index.

        Parameters
        ----------
        context : openmm.Context
        window : int
            Index into self.schedule.
        """
        lambda_coul, lambda_vdw = self.schedule[window]
        context.setParameter('lambda_coul', lambda_coul)
        context.setParameter('lambda_vdw', lambda_vdw)

        for idx in self.mol_atoms:
            q_orig = self._original_charges[idx]
            self._nbforce.setParticleParameters(idx, q_orig * lambda_coul, 1.0, 0.0)
        self._nbforce.updateParametersInContext(context)

    @staticmethod
    def default_lambda_schedule(n_lambda):
        """
        Generate a lambda schedule with n_lambda windows.

        First ~1/3 of windows decouple Coulomb (lambda_coul 1->0, lambda_vdw=1).
        Remaining ~2/3 decouple LJ (lambda_vdw 1->0, lambda_coul=0).

        Parameters
        ----------
        n_lambda : int
            Total number of lambda windows.

        Returns
        -------
        schedule : list of (float, float)
            Each element is (lambda_coul, lambda_vdw).
        """
        n_coul = max(2, int(round(n_lambda / 3)))
        n_vdw = n_lambda - n_coul

        schedule = []
        for i in range(n_coul):
            lam_coul = 1.0 - i / (n_coul - 1) if n_coul > 1 else 0.0
            schedule.append((lam_coul, 1.0))

        for i in range(n_vdw):
            lam_vdw = 1.0 - (i + 1) / n_vdw
            schedule.append((0.0, lam_vdw))

        return schedule

    @staticmethod
    def compute_mbar(u_kn, N_k, temperature):
        """
        Run MBAR analysis on a reduced energy matrix.

        Parameters
        ----------
        u_kn : np.ndarray, shape (n_lambda, n_total)
            Reduced potential energy matrix. u_kn[k, n] = beta * dU to state k
            for sample n.
        N_k : np.ndarray, shape (n_lambda,)
            Number of samples from each state.
        temperature : float
            Temperature in K.

        Returns
        -------
        mu_ex : float
            Excess chemical potential in kJ/mol.
        error : float
            Statistical uncertainty in kJ/mol.
        dG_adj : np.ndarray, shape (n_lambda - 1,)
            Free energy difference between adjacent windows in kJ/mol.
        dG_adj_err : np.ndarray, shape (n_lambda - 1,)
            Uncertainty of dG_adj in kJ/mol.
        """
        import logging
        _logger = logging.getLogger('pymbar')
        _level = _logger.level
        _logger.setLevel(logging.ERROR)
        import pymbar
        _logger.setLevel(_level)

        R = constant.BOLTZMANN * constant.AVOGADRO / 1000  # kJ/(mol*K)
        beta = 1.0 / (R * temperature)

        with np.errstate(under='ignore'):
            mbar = pymbar.MBAR(u_kn, N_k)
            results = mbar.compute_free_energy_differences()

        delta_f = results['Delta_f'][0, -1] / beta
        delta_f_err = results['dDelta_f'][0, -1] / beta

        mu_ex = -delta_f
        error = delta_f_err

        n_lambda = len(N_k)
        dG_adj = np.array([
            results['Delta_f'][i, i + 1] / beta for i in range(n_lambda - 1)
        ])
        dG_adj_err = np.array([
            results['dDelta_f'][i, i + 1] / beta for i in range(n_lambda - 1)
        ])

        return mu_ex, error, dG_adj, dG_adj_err
