import os
import shutil
import math

from .gmx import GmxSimulation
from ...analyzer import is_converged, block_average, average_of_blocks


class NptPPM(GmxSimulation):
    '''
    Npt periodic perturbation simulation with GROMACS enables the calculation of shear viscosity.

    Parameters
    ----------
    packmol : Packmol
    dff : DFF
    gmx : GMX
    jobmanager : subclass of JobManager
    amplitudes_steps : dict, [float, int]

    Attributes
    ----------
    packmol : Packmol
    dff : Packmol
    gmx : GMX
    jobmanager : subclass of JobManager
    n_atom_default : int
    n_mol_default : int
    n_mol_list : list of int
    amplitudes_steps : dict, [float, int]
    logs : list of str
    '''

    def __init__(self, amplitudes_steps=None, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'npt-ppm'
        self.n_atoms_default = 6000
        self.amplitudes_steps = amplitudes_steps or {0.010: int(3.0e6),
                                                     0.020: int(1.5e6),
                                                     0.030: int(1.0e6),
                                                     0.040: int(0.5e6),
                                                     0.050: int(0.5e6)}
        self.logs = ['ppm-%.3f.log' % ppm for ppm in self.amplitudes_steps.keys()]

    def set_system(self, smiles_list, **kwargs):
        super().set_system(smiles_list, **kwargs)

    def build(self, export=True, ppf=None):
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, 'init.pdb', length=self.length - 2,
                               silent=True)

        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd,
                                         mol_corr='init.pdb',
                                         length=self.length)
        if export:
            self.export(ppf=ppf)

    def prepare(self, prior_job_dir=None, gro='conf.gro', top='topol.top', T=298, P=1, jobname=None,
                dt=0.001, nst_eq=int(1E5), nst_edr=50, replicate=None,
                **kwargs) -> [str]:
        if prior_job_dir is None:
            raise Exception('prior_job_dir is needed for PPM simulation')

        # Copy gro and topology files from prior NPT simulation
        shutil.copy(os.path.join(prior_job_dir, gro), '.')
        shutil.copy(os.path.join(prior_job_dir, top), '.')
        for f in os.listdir(prior_job_dir):
            if f.endswith('.itp'):
                shutil.copy(os.path.join(prior_job_dir, f), '.')

        if replicate is not None:
            self.gmx.replicate_gro(gro, top, replicate)

        nprocs = self.jobmanager.nprocs
        commands = []

        for ppm, nst_run in self.amplitudes_steps.items():
            name_eq = 'eq-%.3f' % ppm
            name_ppm = 'ppm-%.3f' % ppm

            # NPT-PPM equilibrium with Nose-Hoover thermostat and Parrinello-Rahman barostat
            # TODO should test the validity of V-rescale thermostat. V-rescale is preferred if it works
            self.gmx.prepare_mdp_from_template('t_npt_ppm.mdp', mdp_out='grompp-%s.mdp' % name_eq,
                                               T=T, P=P,
                                               nsteps=nst_eq, nstxtcout=0, restart=True,
                                               tcoupl='nose-hoover', ppm=ppm)
            cmd = self.gmx.grompp(mdp='grompp-%s.mdp' % name_eq, gro=gro, top=top,
                                  tpr_out='%s.tpr' % name_eq, get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name=name_eq, nprocs=nprocs, get_cmd=True)
            commands.append(cmd)

            # NPT-PPM production with Nose-Hoover thermostat and Parrinello-Rahman barostat
            # TODO should test the validity of V-rescale thermostat. V-rescale is preferred if it works
            self.gmx.prepare_mdp_from_template('t_npt_ppm.mdp', mdp_out='grompp-%s.mdp' % name_ppm,
                                               T=T, P=P,
                                               dt=dt, nsteps=nst_run, nstenergy=nst_edr,
                                               restart=True,
                                               tcoupl='nose-hoover', ppm=ppm)
            cmd = self.gmx.grompp(mdp='grompp-%s.mdp' % name_ppm, gro='%s.gro' % name_eq, top=top,
                                  cpt='%s.cpt' % name_eq, tpr_out='%s.tpr' % name_ppm, get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name=name_ppm, nprocs=nprocs, get_cmd=True)
            commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def run(self):
        super().run()

    def extend(self, extend=500, jobname=None, sh=None) -> [str]:
        '''
        extend simulation for 500 ps
        '''
        nprocs = self.jobmanager.nprocs
        commands = []
        for ppm in self.amplitudes_steps.keys():
            name_ppm = 'ppm-%.3f' % ppm
            self.gmx.extend_tpr(name_ppm, extend, silent=True)
            # Extending NPT production
            cmd = self.gmx.mdrun(name=name_ppm, nprocs=nprocs, extend=True, get_cmd=True)
            commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure, sh=sh)
        return commands

    def analyze(self, dirs=None, check_converge=True):
        import numpy as np
        from ...panedr import edr_to_df
        from ...analyzer.fitting import polyfit

        vis_list = []
        stderr_list = []
        for ppm in self.amplitudes_steps.keys():
            name_ppm = 'ppm-%.3f' % ppm
            df = edr_to_df('%s.edr' % name_ppm)
            inv_series = df['1/Viscosity']

            # select last 4/5 of data
            when = inv_series.index[len(inv_series) // 5]

            # use block average to estimate stderr, because 1/viscosity fluctuate heavily
            inv_blocks = average_of_blocks(inv_series.loc[when:])
            vis_blocks = [1000 / inv for inv in inv_blocks]  # convert Pa*s to cP
            vis_list.append(np.mean(vis_blocks))
            stderr_list.append(np.std(vis_blocks, ddof=1) / math.sqrt(len(vis_blocks)))

        coef_, score = polyfit(self.amplitudes_steps.keys(), vis_list, 1)
        # coef_, score = polyfit(self.amplitudes_steps.keys(), vis_list, 1, weight=1 / np.sqrt(stderr_list))

        return {
            'viscosity'  : coef_[0],
            'score'      : score,
            'vis_list'   : vis_list,
            'stderr_list': stderr_list,
        }

    def clean(self):
        for f in os.listdir(os.getcwd()):
            if f.startswith('eq.') or f.startswith('#'):
                try:
                    os.remove(f)
                except:
                    pass
