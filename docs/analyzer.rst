Analyzer
========

Analyzer module provides a bunch of shortcuts for performing the most commonly used analysing method,
like curve fitting, structure determination and time series analysis.

Structural analysis
-------------------

.. currentmodule:: mstk.analyzer.structure

.. autosummary::
    :toctree: _generated/

    calc_weighted_average
    calc_com
    calc_rg

Structural analysis for vapor-liquid interface
----------------------------------------------

.. currentmodule:: mstk.analyzer.vle

.. autosummary::
    :toctree: _generated/

    check_vle_density
    N_vaporize_condense

Time series analysis
--------------------

.. currentmodule:: mstk.analyzer.series

.. autosummary::
    :toctree: _generated/

    block_average
    average_of_blocks
    is_converged
    efficiency_with_block_size

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
