Installation
============

Install dependencies
--------------------

`mstk` requires `numpy`, `pandas` and `rdkit` for most of its functionalities.
`openmm` is required for energy calculation, minimization etc...
`chemfiles` is required for parsing `XTC` and `DCD` trajectory formats.
`packmol` is required for building simulation box.
`scipy`, `scikit-learn`, `matplotlib-base` and `pymbar-core` are required by `analyzer` module.

It is recommended to install these dependencies with `conda`.

.. code-block:: bash

    conda install -c conda-forge numpy pandas rdkit openmm chemfiles packmol scipy scikit-learn matplotlib-base pymbar-core


Install `mstk`
--------------

`mstk` itself can be installed either with `pip` or `conda`.

Install `mstk` with pip
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    pip install mstk

Install `mstk` with conda
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    conda install -c conda-forge -c z-gong mstk
