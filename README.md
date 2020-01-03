# mstools
A set of tools or wrappers for preparing, running and analyzing molecular simulations  

#### analyzer
Tools for analyzing timeseries. *e.g.* detect equlibriation, detect vapor-liquid separation, perform 1-D or 2-D fitting

#### forcefield
Library for handling force field

#### jobmanager
Wrappers for HPC job manager like `Slurm`, `Torque`

#### omm
Toolkit for simplifying the process of running simulations with OpenMM

#### saved_mol2
Prebuild structures for several weired molecules (PF6 and some amides) which cannot be correcly handled by `OpenBabel`

#### simulation
Predefined simulation protocols for high throughput MD and QM computation
* `gauss/cv.py` for calculating heat capacity using Gaussian
* `gmx/npt.py` for running and analyzing bulk liquid simulation using GROMACS
* `gmx/nvt_slab.py` vapor-liquid coexistence simulation
* `gmx/npt_ppm.py` PPM viscosity simulation
* `gmx/nvt.py` NVT simulation to calculate diffusion and GK viscosity
* `gmx/nvt_gas.py` gas phase in periodic box. Useful for calculating Hvap of ionic liquids
* `gmx/nvt_vacuum.py` gas phase in vacuum
* `gmx/nvt_example.py` basic single molecule example box

#### template
Templates for input files of `GROMACS` and `DFF`

#### trajectory
Library for manipulating and analyzing simulation trajectories

#### wrapper
Wrappers for `GROMACS`, `Packmol`, `DFF`, `Gaussian`  
A parser for manipulate `ppf` file used by DFF is also provided

#### errors, formula, panedr, utils
Misc utils for parsing molecular formula, reading GROMACS energy files and so on

# scripts
Some example scripts for running and analyzing simulations using `mstools`

#### An example of using ms-tools
This script build simulation box of hexane using `Packmol`,
assign temperature-dependent *TEAM_MGI* force field using `DFF`,
prepare input files for running `GROMACS` NPT simulation at 500 K and 10 bar,
and generate script for `Slurm` job manager
```
import sys
sys.path.append(/PATH/OF/MS-TOOLS)

from mstools.wrapper import Packmol, DFF, GMX
from mstools.jobmanager import Slurm
from mstools.simulation.gmx import Npt

packmol = Packmol(packmol_bin='/share/apps/tools/packmol')
dff = DFF(dff_root='/home/gongzheng/apps/DFF/Developing', default_db='TEAMFF', default_table='MGI')
gmx = GMX(gmx_bin='gmx_serial', gmx_mdrun='gmx_gpu mdrun')
slurm = Slurm(*('cpu', 16, 0, 16), env_cmd='module purge; module load gromacs/2016.6')

npt = Npt(packmol=packmol, dff=dff, gmx=gmx, jobmanager=slurm)

npt.set_system(['CCCCCC'])
npt.build()
npt.prepare(T=500, P=10, drde=True)
```
