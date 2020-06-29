
Analyzer
========

Analyzer module provides a bunch of shortcuts for performing the most commonly used analysing method,
like curve fitting, structure determination and time series analysing.

Structural analysis
-------------------

.. currentmodule:: mstools.analyzer.structure

.. autosummary::
    :toctree: _generated/

    check_vle_density
    N_vaporize_condense

Time series analysis
--------------------

.. currentmodule:: mstools.analyzer.series

.. autosummary::
    :toctree: _generated/

    block_average
    average_of_blocks
    is_converged
    efficiency_with_block_size

Curve fitting
-------------

.. currentmodule:: mstools.analyzer.fitting

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
