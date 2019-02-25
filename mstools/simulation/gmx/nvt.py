import os, shutil
from .gmx import GmxSimulation


class Nvt(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'nvt'
        self.requirement = []
        self.logs = []

    def build(self, ppf=None, minimize=False):
        pass

    def prepare(self, prior_job_dir=None, gro='npt.gro', top='topol.top', T=298, jobname=None,
                n_msd=0, nst_msd=int(1E6),
                n_velacc=0, nst_velacc=int(1E5), nstvout_velacc=10,
                n_vis=0, nst_vis=int(5E5), nstenergy_vis=2,
                **kwargs) -> [str]:
        if prior_job_dir is None:
            raise Exception('prior_job_dir is needed for NVT simulation')

        # Copy gro and topology files from prior NPT simulation
        shutil.copy(os.path.join(prior_job_dir, gro), '.')
        shutil.copy(os.path.join(prior_job_dir, top), '.')
        for f in os.listdir(prior_job_dir):
            if f.endswith('.itp'):
                shutil.copy(os.path.join(prior_job_dir, f), '.')
        # Scale gro box for NVT simulation
        box = self.gmx.get_box(os.path.join(prior_job_dir, 'npt.edr'))
        self.gmx.scale_box(gro, gro, box)

        nprocs = self.jobmanager.nprocs
        commands = []

        if n_msd > 0 or n_velacc > 0:
            self.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-eq.mdp', T=T, nsteps=int(5E4),
                                               nstxtcout=0, tcoupl='nose-hoover')
            # equilibrium
            cmd = self.gmx.grompp(mdp='grompp-eq.mdp', gro=gro, top=top, tpr_out='eq.tpr', get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name='eq', nprocs=nprocs, get_cmd=True)
            commands.append(cmd)
        if n_msd > 0:
            self.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-msd.mdp', T=T, nsteps=nst_msd,
                                               nstxtcout=10000, tcoupl='nose-hoover', restart=True)
            cmd = self.gmx.grompp(mdp='grompp-msd.mdp', gro='eq.gro', top=top, tpr_out='nvt-msd.tpr',
                                  cpt='eq.cpt', get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name='nvt-msd', nprocs=nprocs, get_cmd=True)
            commands.append(cmd)
        if n_velacc > 0:
            self.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-velacc.mdp', T=T, nsteps=nst_velacc,
                                               nstvout=nstvout_velacc, tcoupl='nose-hoover', restart=True)
            cmd = self.gmx.grompp(mdp='grompp-velacc.mdp', gro='eq.gro', top=top, tpr_out='nvt-velacc.tpr',
                                  cpt='eq.cpt', get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name='nvt-velacc', nprocs=nprocs, get_cmd=True)
            commands.append(cmd)
        if n_vis > 0:
            self.gmx.prepare_mdp_from_template('t_nvt_anneal.mdp', mdp_out='grompp-anneal.mdp', T=T, nsteps=int(1.5E5))
            self.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-vis.mdp', T=T, nsteps=nst_vis,
                                               nstenergy=nstenergy_vis, tcoupl='nose-hoover', restart=True)
        for i in range(n_vis):
            gro_anneal = 'anneal%i.gro' % i
            name_anneal = 'anneal%i' % i
            name_vis = 'nvt-vis%i' % i
            tpr_anneal = name_anneal + '.tpr'
            tpr_vis = name_vis + '.tpr'

            # equilibrium
            cmd = self.gmx.grompp(mdp='grompp-anneal.mdp', gro=gro, top=top, tpr_out=tpr_anneal, get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name=name_anneal, nprocs=nprocs, get_cmd=True)
            commands.append(cmd)

            # viscosity from Green-Kubo
            cmd = self.gmx.grompp(mdp='grompp-vis.mdp', gro=gro_anneal, top=top, tpr_out=tpr_vis,
                                  cpt=name_anneal + '.cpt', get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name=name_vis, nprocs=nprocs, get_cmd=True)
            commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def analyze(self, dirs=None):
        pass
