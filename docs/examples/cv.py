from mstools.wrapper import Gauss
from mstools.jobmanager import Slurm
from mstools.simulation.gauss import Cv

gauss = Gauss(gauss_bin='/share/apps/g09/g09')
slurm = Slurm(queue='cpu', nprocs=16, ngpu=0, nprocs_request=16,
              env_cmd='module purge; module load gaussian')

cv = Cv(gauss=gauss, jobmanager=slurm)
cv.set_system(['CCCCCC'])
cv.prepare(T_list=[100, 200,300,400,500,600,700], n_conformer=3)
