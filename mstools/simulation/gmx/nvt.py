import os, shutil
import math
from typing import Dict
from .gmx import GmxSimulation
from ...unit import Unit


class Nvt(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'nvt'
        self.requirement = []
        self.logs = []
        for i in range(5):
            self.logs.append('cv%i.log' % i)
            self.logs.append('dos%i.log' % i)

    def build(self, ppf=None, minimize=False):
        pass

    def prepare(self, gro='conf.gro', top='topol.top', T=None, P=None, jobname=None,
                prior_job_dir=None, nst_anneal=int(1E5), nst_cv=int(4E4), nst_vis=int(5E5), n_cv=5, n_vis=25, **kwargs):
        # Copy topology files from prior NPT simulation
        shutil.copy(os.path.join(prior_job_dir, top), '.')
        for f in os.listdir(prior_job_dir):
            if f.endswith('.itp'):
                shutil.copy(os.path.join(prior_job_dir, f), '.')

        # Scale gro box for NVT simulation
        box = self.gmx.get_box(os.path.join(prior_job_dir, 'npt.edr'))
        shutil.copy(os.path.join(prior_job_dir, 'npt.gro'), gro)
        self.gmx.scale_box(gro, gro, box)

        nprocs = self.jobmanager.nprocs
        commands = []

        # Heat capacity using 2-Phase Thermodynamics
        self.gmx.prepare_mdp_from_template('t_nvt_anneal.mdp', mdp_out='grompp-anneal.mdp', T=T, nsteps=nst_anneal)
        self.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-cv.mdp', T=T,
                                           nsteps=nst_cv, nstvout=4, tcoupl='nose-hoover', restart=True)
        self.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-vis.mdp', T=T,
                                           nsteps=nst_vis, nstenergy=2, tcoupl='nose-hoover', restart=True)
        for i in range(max(n_cv, n_vis)):
            gro_anneal = 'anneal%i.gro' % i
            name_anneal = 'anneal%i' % i
            name_cv = 'cv%i' % i
            name_vis = 'vis%i' % i
            tpr_anneal = name_anneal + '.tpr'
            tpr_cv = name_cv + '.tpr'
            tpr_vis = name_vis + '.tpr'

            # equilibrium
            cmd = self.gmx.grompp(mdp='grompp-anneal.mdp', gro=gro, top=top, tpr_out=tpr_anneal, get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name=name_anneal, nprocs=nprocs, get_cmd=True)
            commands.append(cmd)

            # Cv from Density of States
            if i < n_cv:
                cmd = self.gmx.grompp(mdp='grompp-cv.mdp', gro=gro_anneal, top=top, tpr_out=tpr_cv,
                                      cpt=name_anneal + '.cpt', get_cmd=True)
                commands.append(cmd)
                cmd = self.gmx.mdrun(name=name_cv, nprocs=nprocs, get_cmd=True)
                commands.append(cmd)
                cmd = self.gmx.dos(trr=name_cv + '.trr', tpr=tpr_cv, T=T, log_out='dos%i.log' % i, get_cmd=True)
                commands.append(cmd)

            # viscosity from Green-Kubo
            if i < n_vis:
                cmd = self.gmx.grompp(mdp='grompp-vis.mdp', gro=gro_anneal, top=top, tpr_out=tpr_vis,
                                      cpt=name_anneal + '.cpt', get_cmd=True)
                commands.append(cmd)
                cmd = self.gmx.mdrun(name=name_vis, nprocs=nprocs, get_cmd=True)
                commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def analyze(self, dirs=None):
        cv_list = []
        for i in range(5):
            with open('dos%i.log' % i) as f_dos:
                lines = f_dos.readlines()
            for line in lines:
                if line.startswith('Heat capacity'):
                    cv_list.append(float(line.split()[2]))
                    break
            else:
                raise Exception('Heat capacity not found')

        import numpy as np
        return {'cv': [np.mean(cv_list), np.std(cv_list, ddof=1) / math.sqrt(5)]}
