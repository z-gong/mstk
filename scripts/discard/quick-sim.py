#!/usr/bin/env python3

import sys, shutil

from mstools.wrapper import Packmol, DFF, GMX
from mstools.jobmanager import Slurm
from mstools.simulation.gmx import Npt

packmol = Packmol(packmol_bin='/share/apps/tools/packmol')
dff = DFF(dff_root='/home/gongzheng/apps/DFF/Developing', default_db='TEAMFF', default_table='MGI')
gmx = GMX(gmx_bin='gmx_serial', gmx_mdrun='gmx_gpu mdrun')
slurm = Slurm(*('cpu', 16, 0, 16), env_cmd='module purge; module load gromacs/2016.6')

npt = Npt(packmol=packmol, dff=dff, gmx=gmx, jobmanager=slurm)

smiles = sys.argv[1]
T = int(sys.argv[2])
P = int(sys.argv[3])
jobname = '%s-%i-%i' % (smiles, T, P)

if len(sys.argv) > 4:
    ppf = 'ff.ppf'
    shutil.copy(sys.argv[4], ppf)
else:
    ppf = None

npt.set_system([smiles], n_atoms=3000)
npt.build(export=True, ppf=ppf)
npt.prepare(drde=True, T=T, P=P, TANNEAL=None, nst_eq=int(2E5), nst_run=int(3E5), nst_xtc=500, jobname=jobname)
