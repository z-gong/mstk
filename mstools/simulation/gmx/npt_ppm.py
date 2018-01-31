import os
import shutil

from .gmx import GmxSimulation
from ...analyzer import is_converged, ave_and_stderr


class NptPPM(GmxSimulation):
    def __init__(self, amplitudes=(0.02, 0.04, 0.06, 0.08, 0.10), **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'nvt-ppm'
        self.requirement = []
        self.logs = []
        self.amplitudes = amplitudes

    def build(self, export=True, ppf=None, minimize=False):
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, 'init.pdb', length=self.length - 2, silent=True)

        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd, mol_corr='init.pdb',
                                         length=self.length)
        if export:
            self.export(ppf=ppf, minimize=minimize)

    def prepare(self, prior_job_dir=None, gro='conf.gro', top='topol.top', T=298, P=1, jobname=None,
                dt=0.002, nst_eq=int(1E5), nst_run=int(1E6), nst_edr=50,
                **kwargs) -> [str]:
        if prior_job_dir == None:
            raise Exception('prior_job_dir is needed for PPM simulation')

        # Copy gro and topology files from prior NPT simulation
        shutil.copy(os.path.join(prior_job_dir, gro), '.')
        shutil.copy(os.path.join(prior_job_dir, top), '.')
        for f in os.listdir(prior_job_dir):
            if f.endswith('.itp'):
                shutil.copy(os.path.join(prior_job_dir, f), '.')

        nprocs = self.jobmanager.nprocs
        commands = []

        for ppm in self.amplitudes:
            name_eq = 'eq-%.2f' % ppm
            name_ppm = 'ppm-%.2f' % ppm

            # NPT-PPM equilibrium with Nose-Hoover thermostat and Berendsen barostat
            self.gmx.prepare_mdp_from_template('t_npt_ppm.mdp', mdp_out='grompp-%s.mdp' % name_eq, T=T, P=P,
                                               nsteps=nst_eq, nstxtcout=0, restart=True,
                                               tcoupl='nose-hoover', ppm=ppm)
            cmd = self.gmx.grompp(mdp='grompp-%s.mdp' % name_eq, gro=gro, top=top,
                                  tpr_out='%s.tpr' % name_eq, get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name=name_eq, nprocs=nprocs, get_cmd=True)
            commands.append(cmd)

            # NPT-PPM production with Nose-Hoover thermostat and Parrinello-Rahman barostat
            self.gmx.prepare_mdp_from_template('t_npt_ppm.mdp', mdp_out='grompp-%s.mdp' % name_ppm, T=T, P=P,
                                               dt=dt, nsteps=nst_run, nstenergy=nst_edr, restart=True,
                                               tcoupl='nose-hoover', ppm=ppm)
            cmd = self.gmx.grompp(mdp='grompp-%s.mdp' % name_ppm, gro='%s.gro' % name_eq, top=top,
                                  cpt='%s.cpt' % name_eq, tpr_out='%s.tpr' % name_ppm, get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name=name_ppm, nprocs=nprocs, get_cmd=True)
            commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def extend(self, extend=500, jobname=None, sh=None) -> [str]:
        '''
        extend simulation for 500 ps
        '''
        nprocs = self.jobmanager.nprocs
        commands = []
        for ppm in self.amplitudes:
            name_ppm = 'ppm-%.2f' % ppm
            self.gmx.extend_tpr(name_ppm, extend)
            # Extending NPT production
            cmd = self.gmx.mdrun(name=name_ppm, nprocs=nprocs, extend=True, get_cmd=True)
            commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure, sh=sh)
        return commands

    def analyze(self, dirs=None):
        from ...panedr import edr_to_df
        from ...analyzer.fitting import polyfit
        import numpy as np

        vis_list = []
        stderr_list = []
        for ppm in self.amplitudes:
            name_ppm = 'ppm-%.2f' % ppm
            df = edr_to_df('%s.edr' % name_ppm)
            density_series = df.Density
            inv_vis_series = df['1/Viscosity']

            converged, when = is_converged(density_series)

            if not converged:
                when = density_series.index[len(density_series) // 2]

            inv_vis, stderr = ave_and_stderr(inv_vis_series.loc[when:])
            vis_list.append(1000 / inv_vis)  # convert Pa.s to cP
            stderr_list.append(stderr)

        # coef_ = np.polynomial.polynomial.polyfit(self.amplitudes, vis_list, 1, w=stderr_list)
        coef_, score = polyfit(self.amplitudes, vis_list, 1)

        return {
            'viscosity': coef_[0],
            'score': score,
            'vis_list': vis_list,
        }

    def clean(self):
        for f in os.listdir(os.getcwd()):
            if f.startswith('eq.') or f.startswith('#'):
                try:
                    os.remove(f)
                except:
                    pass
