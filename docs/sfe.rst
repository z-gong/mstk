
Solvation Free Energy
=====================

Alchemical coupling for computing solvation free energy (excess chemical potential)
using OpenMM. The molecule is coupled into its environment in two phases:
LJ interactions are turned on first using a soft-core potential, then Coulomb
interactions are turned on (with LJ fully on).

Lambda Schedule
---------------

- **Phase 1** (~3/4 of windows): ``lambda_vdw`` 0→1 (soft-core), ``lambda_coul`` = 0.
- **Phase 2** (~1/4 of windows): ``lambda_coul`` 0→1, ``lambda_vdw`` = 1.

Initial Force Decomposition
----------------------------

``system.to_omm_system()`` produces (bonded forces omitted):

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - Force
     - Computes
     - Per-particle params
   * - NonbondedForce (PME)
     - Coulomb (real + reciprocal + self)
     - charge, sigma=1.0, eps=0
   * - NonbondedForce exceptions
     - 1-4 Coulomb + 1-2/1-3 exclusions
     - chargeProd
   * - CustomNonbondedForce (LJ)
     - ``A(type1,type2)*invR6^2 - B(type1,type2)*invR6``
     - type (int)
   * - CustomBondForce (1-4 LJ)
     - ``C*eps*((sigma/r)^n - (sigma/r)^m)``
     - eps, sigma, n, m

Final Alchemical Force Table
----------------------------

After modification (M = decoupled molecule, env = everything else):

.. list-table::
   :header-rows: 1
   :widths: 30 30 20

   * - Force
     - Scope
     - Lambda-dependent?
   * - NonbondedForce (PME)
     - All pairs (unchanged)
     - M charges × lambda_coul
   * - NonbondedForce exceptions
     - 1-4 Coulomb + exclusions (unchanged)
     - No
   * - Original CustomNonbondedForce
     - (env, env) + (M, M)
     - No
   * - Alchemical CustomNonbondedForce
     - (M, env) soft-core
     - lambda_vdw
   * - Correction CustomBondForce
     - M-M beyond-1-4 pairs
     - lambda_coul
   * - 1-4 CustomBondForce
     - All 1-4 LJ (unchanged)
     - No

API Reference
-------------

.. currentmodule:: mstk.sfe

.. autosummary::
    :toctree: _generated/

    SFEManager
