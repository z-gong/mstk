
Trajectory
==========

.. currentmodule:: mstools.trajectory

.. autosummary::
    :toctree: _generated/

    Trajectory
    Frame

Trajectory handler
------------------

Handlers for reading/writing different trajectory formats by :class:`~mstools.trajectory.Trajectory`.
They are not meant to be called directly.
However, you can pass a handler class when initializing a :class:`~mstools.trajectory.Trajectory` object.

.. currentmodule:: mstools.trajectory

.. autosummary::
    :toctree: _generated/

    TrjHandler
    Dcd
    Gro
    LammpsTrj
    Xtc
    Xyz
    CombinedTrj
