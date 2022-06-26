from mstk.wrapper import Gauss
from mstk.jobmanager import Slurm
from mstk.simulation.gauss import Cv

# Initialize wrappers for Gaussian 09
gauss = Gauss(gauss_bin='/share/apps/g09/g09')

# Initialize Slurm job scheduler.
# This calculation will be submitted to 'cpu' partition.
# 16 CPU cores without hyper-threading will be consumed.
# The module environment will be configured before running Gaussian.
slurm = Slurm(queue='cpu', nprocs=16, ngpu=0, nprocs_request=16,
              env_cmd='module purge; module load gaussian')

# Initialize a Cv QM calculation on dodecane.
cv = Cv(gauss=gauss, jobmanager=slurm)
cv.set_system(['CCCCCCCCCCCC'])

# Generate input files for running Cv QM calculation.
# Three random conformers will be created for calculation.
# The thermodynamic analysis will be performed at several temperature points.
cv.prepare(T_list=[100, 200,300,400,500,600,700], n_conformer=3)
