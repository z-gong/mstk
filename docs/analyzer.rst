Analyzer
========

Analyzer module provides a bunch of shortcuts for performing the most commonly used analysing method,
like curve fitting and structure determination.

.. note::
   *Deprecated. This module is developed for personal projects. It's not a comprehensive analysis toolkit.*

Structural analysis
-------------------

.. currentmodule:: mstk.analyzer.structure

.. autosummary::
    :toctree: _generated/

    calc_weighted_average
    calc_com
    calc_rg
    calc_hull_volume

Structural analysis for vapor-liquid interface
----------------------------------------------

.. currentmodule:: mstk.analyzer.vle

.. autosummary::
    :toctree: _generated/

    check_vle_density
    N_vaporize_condense

Curve fitting
-------------

.. currentmodule:: mstk.analyzer.fitting

.. autosummary::
    :toctree: _generated/

    polyfit
    polyval
    polyval_derivative
    polyfit_2d
    polyval_derivative_2d
    curve_fit_rsq
    fit_vle_dminus
    fit_vle_dplus
