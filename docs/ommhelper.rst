
OpenMM helper
=============

This module provides various tools to simplify the usage or extend the functionality of OpenMM.
Mainly, they can be categorized into file parsers, pre-constructed forces and reporters.

Parsers
-------

.. currentmodule:: mstools.ommhelper

.. autosummary::
    :toctree: _generated/

    OplsPsfFile
    GroFile


Forces
------

.. currentmodule:: mstools.ommhelper.force

.. autosummary::
    :toctree: _generated

    slab_correction
    restrain_particle_number
    wall_lj126
    wall_power
    electric_field
    spring_self
    CLPolCoulTT

Reporters
---------

.. currentmodule:: mstools.ommhelper.reporter

.. autosummary::
    :toctree: _generated/

    ViscosityReporter
    DrudeTemperatureReporter
    GroReporter
    CheckpointReporter
    StateDataReporter
