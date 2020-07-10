Examples
========

Assign atom types using typing engine
-------------------------------------

`mstools` provide two typing engines for assigning atom types.
The first one :class:`~mstools.forcefield.typer.ZftTyper` is based on local chemical environment
described by SMARTS pattern and a hierarchical rule.
A type definition should be fed to this typing engine, either by file name or by string.
In order to construct the type definition file correctly, knowledge about SMARTS is required.
https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html.
The idea of the hierarchical rule is described in http://dx.doi.org/10.1002/jcc.24244.
Here is an example of assigning OPLS atom types for butane and benzene.

.. literalinclude:: examples/type.py

Export input files for simulation engines
-----------------------------------------

Use pre-defined simulation protocol
-----------------------------------

`mstools` provides several pre-defined protocols for fast modeling of small molecules.
It is initially designed for high-throughput prediction of liquids properties with TEAM_MGI force field.
Therefore, the requirement for some packages are hardcoded.
Packmol, DFF and GROMACS should have been installed.

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
* _job_slurm.sh


This is an example to use the pre-defined Gaussian protocol to calculate
the intra-molecular heat capacitance of dodecane.
Gaussian should have been installed.
Considering that dodecane is a flexible molecule, three random conformers are generated.
Though usually conformations have little effect on intra-molecular vibration modes of weak polar molecules.

.. literalinclude:: examples/cv.py

This script will generate four files in current directory:

* conf-0.gjf
* conf-1.gjf
* conf-2.gjf
* _job_slurm.sh

The three gjf files are the input for Gaussian of three random dodecane conformers.
And _job_slurm.sh is the control script for Slurm job scheduler.
