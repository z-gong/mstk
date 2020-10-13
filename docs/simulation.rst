
Simulation protocols
====================

This module provides predefined simulation protocols for performing different types of MD and QM calculations.
e.g. NPT simulation or periodic perturbation simulation with GROMACS, normal mode analysis with Gaussian, etc...
They are designed for high-throughput prediction of liquid properties.

GROMACS protocols
-----------------

Currently, these MD simulation protocols are hardcoded with DFF typing engine and TEAM force field.
It will be subject to refactor in the future to support different typing engines and force fields.

.. currentmodule:: mstools.simulation.gmx

.. autosummary::
    :toctree: _generated/

    GmxSimulation
    Npt
    NvtSlab
    NptPPM
    Nvt
    NvtVacuum
    NvtGas

Gaussian protocols
------------------

.. currentmodule:: mstools.simulation.gauss

.. autosummary::
    :toctree: _generated/

    GaussSimulation
    Cv
