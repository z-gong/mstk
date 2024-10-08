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

`mstk` itself can be installed either with `pip` or from source code.

Install `mstk` with pip
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    pip install mstk

Install `mstk` from source code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clone it from GitHub:

.. code-block:: bash

    git clone https://github.com/z-gong/mstk

Add the following line to .bashrc so that it can be found by python interpreter:

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:/path/of/mstk
    export PATH=$PATH:/path/of/mstk/scripts
