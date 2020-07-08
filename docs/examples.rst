Examples
========

GROMACS protocol
----------------

This is an example to use the pre-defined GROMACS protocol to perform NPT
bulk liquid simulation on hexane.

.. literalinclude:: examples/npt.py

This script will generate following files in current directory:

* conf.gro
* topol.itp
* topol.top
* topol-hvap.top
* em.mdp
* anneal.mdp
* eq.mdp
* npt.mdp
* _job.slurm

Gaussian protocol
-----------------

This is an example to use the pre-defined Gaussian protocol to calculate
the intra-molecular heat capacitance of dodecane.
Considering that dodecane is a flexible molecule, three random conformers are generated.
Though usually conformations have little effect on intra-molecular vibration modes of weak polar molecules.

.. literalinclude:: examples/cv.py

This script will generate four files in current directory:

* conf-0.gjf
* conf-1.gjf
* conf-2.gjf
* _job.slurm

The three gjf files are the input for Gaussian of three random dodecane conformers.
And _job.slurm is the control script for Slurm job scheduler.
