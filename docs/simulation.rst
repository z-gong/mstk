
Simulation protocols
====================

Predefined simulation protocols for performing different types of MD and QM calculations.
They are designed for High-throughput prediction of liquid properties with TEAM force field.
The force field typing and parameter assigning are hardcoded with DFF,
therefore this module is not supposed to be used for other purposes.

It may be subject to refactor in the future to support different typing and simulation engines.

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
    NvtVacuum

Gaussian protocols
------------------

.. currentmodule:: mstools.simulation.gauss

.. autosummary::
    :toctree: _generated/

    GaussSimulation
    Cv
