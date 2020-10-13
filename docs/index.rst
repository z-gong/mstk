.. mstools documentation master file, created by
   sphinx-quickstart on Mon Jun 29 12:36:37 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

mstools
=======

`mstools` is a toolkit for preparing, running and analyzing molecular simulations.
The main functionality of `mstools` can be separated into two parts:

1. The first is for general applications. It simplifies the workflow of setting up a simulation by
manageable atom type and force field assignment, topology manipulation, trajectory reading/writing, etc...

2. The second is for high-throughput molecular simulations. It is designed for predicting liquid properties in a reproducible way.
Pre-defined simulation protocols enables fast and reliable simulation of target molecules by automatic force field assignment,
input generation, error detection, result retrieving and post-processing.

.. toctree::
   :maxdepth: 2

   installation
   examples

Modules
-------

.. toctree::
   :maxdepth: 2

   topology
   forcefield
   simsys
   jobmanager
   simulation
   ommhelper
   trajectory
   analyzer
   wrapper
   misc


License
-------

mstools is licensed under LGPL v2.1