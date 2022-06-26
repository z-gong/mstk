
Simulation system
=================

A simulation system is a combination of topology and force field.
It defines a system that is ready for simulation.

After a system is constructed, it can be exported to different simulation engines for calculation.
Current, LAMMPS, GROMACS, NAMD and OpenMM are supported.
Depending on the functional form used in the force field, some simulation engines may not be supported for specific system.

.. currentmodule:: mstk.simsys

.. autosummary::
    :toctree: _generated/

    System

Exporter
--------

.. currentmodule:: mstk.simsys

.. autosummary::
    :toctree: _generated/

    GromacsExporter
    LammpsExporter
    NamdExporter
    OpenMMExporter
