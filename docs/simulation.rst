
Simulation protocols
====================

This module provides predefined simulation protocols for performing different types of MD and QM calculations.
e.g. NPT simulation or periodic perturbation simulation with GROMACS, normal mode analysis with Gaussian, etc...
They are designed for high-throughput prediction of liquid properties.

Currently, these simulation protocols are hardcoded with DFF typing engine and TEAM force field.
It will be subject to refactor in the future to support different typing engines and force fields.

.. currentmodule:: mstools.simulation

.. autosummary::
    :toctree: _generated/

    Simulation

GROMACS protocols
-----------------

.. currentmodule:: mstools.simulation.gmx

.. autosummary::
    :toctree: _generated/

    GmxSimulation
    Npt
    NptPPM
    Nvt
    NvtSlab
    NvtGas

Gaussian protocols
------------------

.. currentmodule:: mstools.simulation.gauss

.. autosummary::
    :toctree: _generated/

    GaussSimulation
    Cv
