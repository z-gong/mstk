# mstools

A toolkit for preparing, running and analyzing molecular simulations.

## Capabilities
   
  * Assign atom types and force field parameters based on local chemical environment
  * Generate input files for simulation engines like LAMMPS, GROMACS, OpenMM
  * Predict liquid properties in an efficient and reproducible way through pre-defined simulation protocols
  * Access local and remote job schedulers for high-throughput simulations
  * Read/write common topology files (LAMMPS, PSF, ZMAT etc...)
  * Read/write common trajectory files (LAMMPS, GRO, DCD, XTC, etc...)
  

## Example

This script will build a simulation box of hexane using `Packmol`,
assign temperature-dependent *TEAM_MGI* force field using `DFF`,
prepare input files for running `GROMACS` NPT simulation at 500 K and 10 bar,
and generate corresponding script for `Slurm` job scheduler.
```
from mstools.wrapper import Packmol, DFF, GMX
from mstools.jobmanager import Slurm
from mstools.simulation.gmx import Npt

packmol = Packmol(packmol_bin='/share/apps/tools/packmol')
dff = DFF(dff_root='/home/gongzheng/apps/DFF/Developing', default_db='TEAMFF', default_table='MGI')
gmx = GMX(gmx_bin='gmx_serial', gmx_mdrun='gmx_gpu mdrun')
slurm = Slurm(queue='cpu', nprocs=16, ngpu=0, nprocs_request=16,
              env_cmd='module purge; module load gromacs/2016.6')

npt = Npt(packmol=packmol, dff=dff, gmx=gmx, jobmanager=slurm)
npt.set_system(['CCCCCC'])
npt.build()
npt.prepare(T=500, P=10, drde=True)
```

## Documentation

https://ms-tools.readthedocs.io/en/latest/index.html
