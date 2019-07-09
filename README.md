# ms-tools
A set of tools or wrappers for preparing, running and analyzing molecular simulations  

### analyzer
Tools for analyzing timeseries. e.g. detect equlibriation, detect vapor-liquid separation, perform 1-D or 2-D fitting

### jobmanager
Wrappers for HPC job manager like `Slurm`, `Torque`

### wrapper
Wrappers for `GROMACS`, `Packmol`, `DFF`, `Gaussian`  
A parser for manipulate `ppf` file used by DFF is also provided

### simulation
Predefined simulation protocols for high throughput MD and QM computation
* `gauss/cv.py` for calculating heat capacity using Gaussian
* `gmx/npt.py` for running and analyzing bulk liquid simulation using GROMACS
* `gmx/nvt_slab.py` vapor-liquid coexistence simulation
* `gmx/npt_ppm.py` PPM viscosity simulation
* `gmx/nvt.py` NVT simulation to calculate diffusion and GK viscosity
* `gmx/nvt_gas.py` gas phase in periodic box. Useful for calculating Hvap of ionic liquids
* `gmx/nvt_vacuum.py` gas phase in vacuum
* `gmx/nvt_example.py` basic single molecule example box

### template
Templates for input files of `GROMACS` and `DFF`

### saved_mol2
Prebuild structures for several weired molecules (PF6 and some amides)

### errors, formula, panedr, utils
Misc utils for parsing molecular formula, reading GROMACS energy files and so on

## An example of using ms-tools
This script build NPT simulation box of hexane at 500 K and 10 bar using `Packmol` and `DFF` and prepare input files for `GROMACS` and script for `Slurm` job manager
```
import sys
sys.path.append(/PATH/OF/AIMS_Tools)

from mstools.simulation.gmx import Npt
from mstools.jobmanager import Slurm

slurm=Slurm(*('cpu', 16, 0, 16), env_cmd='module purge; module load gromacs/2016.6')
npt = Npt(packmol_bin='/share/apps/tools/packmol',
          dff_root='/home/gongzheng/apps/DFF/Developing',
          dff_db='TEAMFF',
          dff_table='MGI',
          gmx_bin='gmx_serial',
          gmx_mdrun='gmx_gpu mdrun',
          jobmanager=slurm)

smiles = 'CCCCCC'
T = 500
P = 10
jobname='%s-%i-%i' %(smiles, T, P)

npt.set_system([smiles])
npt.build()
npt.prepare(drde=True, T=T, P=P, nst_eq=int(2E5), nst_run=int(3E5), nst_xtc=500, jobname=jobname)
```
