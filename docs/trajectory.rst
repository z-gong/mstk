
Trajectory
==========

.. currentmodule:: mstk.trajectory

.. autosummary::
    :toctree: _generated/

    Trajectory
    Frame

Trajectory handler
------------------

Handlers for reading/writing different trajectory formats by :class:`~mstk.trajectory.Trajectory`.
They are not meant to be called directly.
However, you can pass a handler class when initializing a :class:`~mstk.trajectory.Trajectory` object.

.. currentmodule:: mstk.trajectory

.. autosummary::
    :toctree: _generated/

    TrjHandler
    Dcd
    Gro
    LammpsTrj
    Xtc
    Xyz
    CombinedTrj
