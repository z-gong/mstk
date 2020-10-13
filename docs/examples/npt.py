from mstools.wrapper import Packmol, DFF, GMX
from mstools.jobmanager import Slurm
from mstools.simulation.gmx import Npt

# Initialize wrappers for Packmol, DFF and GROMACS
packmol = Packmol(packmol_bin='/share/apps/tools/packmol')
dff = DFF(dff_root='/share/apps/DFF/Developing', default_db='TEAMFF', default_table='MGI')
gmx = GMX(gmx_bin='gmx_serial', gmx_mdrun='gmx_gpu mdrun')

# Initialize Slurm job scheduler.
# This simulation will be submitted to 'gtx' partition.
# 1 GPU card and 8 CPU cores without hyper-threading will be consumed.
# The module environment will be configured before running GROMACS.
slurm = Slurm(queue='gtx', nprocs=8, ngpu=1, nprocs_request=8,
              env_cmd='module purge; module load gromacs/2016.6')

# Initialize a Npt simulation, which contains pure hexane.
# There will be approximately 3000 atoms in the simulation box.
npt = Npt(packmol=packmol, dff=dff, gmx=gmx, jobmanager=slurm)
npt.set_system(['CCCCCC'], n_atoms=3000)

# Build the initial configuration for simulation.
npt.build()

# Generate input files for running Npt simulation at 500 K and 10 bar.
# The temperature-dependence of LJ parameters will be applied.
npt.prepare(T=500, P=10, drde=True)
