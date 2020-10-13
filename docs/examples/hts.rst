High-throughput simulations
===========================

`mstools` provides several pre-defined protocols for efficient and reproducible modeling of small molecules.
It is designed for high-throughput prediction of liquids properties with `TEAM` force field.
Therefore, the requirement for some programs are hardcoded for now.

Use pre-defined Gromacs protocols
---------------------------------

Note: This example requires that `Packmol` and `DFF` are installed.

Below is an example of using the pre-defined GROMACS `Npt` protocol to
simulate bulk liquid of hexane with the temperature-dependent `TEAM_MGI` force field.
You may refer to :class:`~mstools.simulation.gmx.Npt` documentation for details of options.

.. literalinclude:: npt.py


This script will generate following files in current directory:

* conf.gro (initial configuration)
* topol.itp (force field parameters)
* topol.top (topology file)
* topol-hvap.top (topology for calculating cohesive energy)
* grompp-em.mdp (control script for energy minimization)
* grompp-anneal.mdp (control script for annealing)
* grompp-eq.mdp (control script for initial equilibration)
* grompp-npt.mdp (control script for production run)
* _job_slurm.sh (input script for Slurm job scheduler)

The script `_job_slurm.sh` can be submitted to Slurm job scheduler to start the simulation.

Use pre-defined Gaussian protocols
----------------------------------

This is an example of using the pre-defined Gaussian `Cv` protocol to calculate
the intra-molecular heat capacitance of dodecane.

Considering that dodecane is a flexible molecule, three random conformers are generated.
Though usually conformations have little effect on intra-molecular vibration modes of weak polar molecules.
You may refer to :class:`~mstools.simulation.gauss.Cv` documentation for details of options.

.. literalinclude:: cv.py

This script will generate following files in current directory:

* conf-0.gjf
* conf-1.gjf
* conf-2.gjf
* _job_slurm.sh

The three `GJF` files are the input for Gaussian of three dodecane conformers.
The script `_job_slurm.sh` can be submitted to Slurm job scheduler to start the calculation.
