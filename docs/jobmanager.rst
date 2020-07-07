
Job scheduler
=============

Molecular simulation is computationally intensive.
Therefore, it is usually impractical to run simulations on local computer.
Instead, simulation jobs are usually submitted to a job scheduler on a distributed computing infrastructure.
Then the jobs will be spread and executed on the compute nodes.

Wrappers for several job schedulers are provided in `mstools`.
They are designed to simplify the most common operations of job schedulers,
like querying, submitting and deleting jobs, generate job scripts, etc...
It also provide a cache to the scheduler so that when high-throughput simulations are performed,
the job scheduler will not experience much pressure.

.. currentmodule:: mstools.jobmanager

.. autosummary::
    :toctree: _generated/

    PbsJob
    JobManager
    Slurm
    RemoteSlurm
    Torque
