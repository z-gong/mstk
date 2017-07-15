import os
import shutil

from .gmx import GmxSimulation
from ...analyzer import is_converged, block_average
from ...unit import Unit


class Npt(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'npt'
        self.requirement = []
        self.logs = ['npt.log', 'hvap.log']

    def build(self, ppf=None, minimize=False):
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, 'init.pdb', length=self.length - 2, silent=True)

        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd, mol_corr='init.pdb',
                                         length=self.length)
        self.export(ppf=ppf, minimize=minimize)

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=None, P=None, jobname=None,
                dt=0.001, nst_eq=int(4E5), nst_run=int(5E5), nst_edr=100, nst_trr=int(1E4), nst_xtc=int(1E3), **kwargs):
        if os.path.abspath(model_dir) != os.getcwd():
            shutil.copy(os.path.join(model_dir, gro), gro)
            shutil.copy(os.path.join(model_dir, top), top)
            for f in os.listdir(model_dir):
                if f.endswith('.itp'):
                    shutil.copy(os.path.join(model_dir, f), '.')

        nprocs = self.jobmanager.nprocs
        commands = []
        # energy minimization
        self.gmx.prepare_mdp_from_template('t_em.mdp', mdp_out='grompp-em.mdp')
        cmd = self.gmx.grompp(mdp='grompp-em.mdp', gro=gro, top=top, tpr_out='em.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='em', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NVT annealing from 0 to 2000 K to target T with Langevin thermostat
        self.gmx.prepare_mdp_from_template('t_nvt_anneal.mdp', mdp_out='grompp-anneal.mdp', T=T,
                                           nsteps=int(1E5), nstxtcout=0)
        cmd = self.gmx.grompp(mdp='grompp-anneal.mdp', gro='em.gro', top=top, tpr_out='anneal.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='anneal', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NPT equilibrium with Langevin thermostat and Berendsen barostat
        self.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-eq.mdp', T=T, P=P / Unit.bar,
                                           dt=0.001, nsteps=nst_eq, nstxtcout=0, restart=True, pcoupl='berendsen')
        cmd = self.gmx.grompp(mdp='grompp-eq.mdp', gro='anneal.gro', top=top, tpr_out='eq.tpr',
                              cpt='anneal.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='eq', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NPT production with Langevin thermostat and Parrinello-Rahman barostat
        self.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-npt.mdp', T=T, P=P / Unit.bar,
                                           dt=dt, nsteps=nst_run, nstenergy=nst_edr, nstxout=nst_trr, nstvout=nst_trr,
                                           nstxtcout=nst_xtc, restart=True)
        cmd = self.gmx.grompp(mdp='grompp-npt.mdp', gro='eq.gro', top=top, tpr_out='npt.tpr',
                              cpt='eq.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='npt', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # Rerun enthalpy of vaporization
        commands.append('export GMX_MAXCONSTRWARN=-1')

        top_hvap = 'topol-hvap.top'
        self.gmx.generate_top_for_hvap(top, top_hvap)
        self.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-hvap.mdp', nstxtcout=0, restart=True)
        cmd = self.gmx.grompp(mdp='grompp-hvap.mdp', gro='eq.gro', top=top_hvap, tpr_out='hvap.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='hvap', nprocs=nprocs, rerun='npt.xtc', get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def extend(self, extend=500, jobname=None):
        '''
        extend simulation for 500 ps
        '''
        self.gmx.extend_tpr('npt.tpr', extend)

        nprocs = self.jobmanager.nprocs
        commands = []
        # Extending NPT production with Velocity Rescaling thermostat and Parrinello-Rahman barostat
        cmd = self.gmx.mdrun(name='npt', nprocs=nprocs, extend=True, get_cmd=True)
        commands.append(cmd)

        # Rerun enthalpy of vaporization
        cmd = self.gmx.mdrun(name='hvap', nprocs=nprocs, rerun='npt.xtc', get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)

    def analyze(self, dirs=None):
        if dirs is None:
            dirs = ['.']
        import pandas as pd
        import panedr

        # TODO check structure freezing
        temp_series = pd.Series()
        press_series = pd.Series()
        potential_series = pd.Series()
        density_series = pd.Series()
        e_inter_series = pd.Series()
        for dir in dirs:
            df = panedr.edr_to_df(os.path.join(dir, 'npt.edr'))
            temp_series = temp_series.append(df.Temperature)
            press_series = press_series.append(df.Pressure)
            potential_series = potential_series.append(df.Potential)
            density_series = density_series.append(df.Density)

            df = panedr.edr_to_df(os.path.join(dir, 'hvap.edr'))
            e_inter_series = e_inter_series.append(df.Potential)

        converged, when = is_converged(density_series)

        if converged:
            return {
                'simulation_length': density_series.index[-1],
                'converged_from': when,
                'temperature': block_average(temp_series.loc[when:]),
                'pressure': block_average(press_series.loc[when:]),
                'potential': block_average(potential_series.loc[when:]),
                'density': block_average(density_series.loc[when:]),
                'e_inter': block_average(e_inter_series.loc[when:]),
            }
        else:
            return None
