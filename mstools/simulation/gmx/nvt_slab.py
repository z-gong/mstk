import os
import shutil

from .gmx import GmxSimulation
from ...wrapper.ppf import delta_ppf


class NvtSlab(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'nvt-slab'
        self.requirement = []
        self.logs = ['nvt.log']

    def build(self, export=True, ppf=None, minimize=False):
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, 'init.pdb',
                               size=[self.length - 2, self.length - 2, self.length * 5 - 2],
                               slab=True, silent=True)

        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd, mol_corr='init.pdb',
                                         size=[self.length, self.length, self.length * 5])
        if export:
            self.export(ppf=ppf, minimize=minimize)

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=298, jobname=None,
                dt=0.002, nst_eq=int(4E5), nst_run=int(5E6), nst_edr=100, nst_trr=int(5E4), nst_xtc=int(1E3),
                drde=False, **kwargs) -> [str]:
        if not drde:
            if os.path.abspath(model_dir) != os.getcwd():
                shutil.copy(os.path.join(model_dir, gro), gro)
                shutil.copy(os.path.join(model_dir, top), top)
                for f in os.listdir(model_dir):
                    if f.endswith('.itp'):
                        shutil.copy(os.path.join(model_dir, f), '.')
        else:
            ### Temperature dependent parameters
            if os.path.abspath(model_dir) != os.getcwd():
                shutil.copy(os.path.join(model_dir, self.msd), self.msd)
            delta_ppf(os.path.join(model_dir, 'ff.ppf'), 'ff.ppf', T)
            self.export(ppf='ff.ppf')

        nprocs = self.jobmanager.nprocs
        commands = []

        # energy minimization
        self.gmx.prepare_mdp_from_template('t_em.mdp', mdp_out='grompp-em.mdp')
        cmd = self.gmx.grompp(mdp='grompp-em.mdp', gro=gro, top=top, tpr_out='em.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='em', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NVT equilibrium
        self.gmx.prepare_mdp_from_template('t_nvt_slab.mdp', mdp_out='grompp-eq.mdp', T=T,
                                           nsteps=nst_eq, nstxtcout=0)
        cmd = self.gmx.grompp(mdp='grompp-eq.mdp', gro='em.gro', top=top, tpr_out='eq.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='eq', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NVT production
        self.gmx.prepare_mdp_from_template('t_nvt_slab.mdp', mdp_out='grompp-nvt.mdp', T=T,
                                           dt=dt, nsteps=nst_run, nstenergy=nst_edr, nstxout=nst_trr, nstvout=nst_trr,
                                           nstxtcout=nst_xtc, restart=True)
        cmd = self.gmx.grompp(mdp='grompp-nvt.mdp', gro='eq.gro', top=top, tpr_out='nvt.tpr',
                              cpt='eq.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='nvt', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def extend(self, extend=5000, jobname=None, sh=None) -> [str]:
        '''
        extend simulation for 5000 ps
        '''
        self.gmx.extend_tpr('nvt.tpr', extend)

        nprocs = self.jobmanager.nprocs
        commands = []
        # Extending NVT production
        cmd = self.gmx.mdrun(name='nvt', nprocs=nprocs, extend=True, get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure, sh=sh)
        return commands

    def analyze(self, dirs=None):
        import pandas as pd
        from ...analyzer.series import is_converged
        from ...analyzer.structure import is_interface
        from ...panedr import edr_to_df

        if dirs is None:
            dirs = ['.']

        df = edr_to_df('nvt.edr')
        potential_series = df.Potential

        converged, when = is_converged(potential_series, frac_min=0.25)
        if not converged:
            return None

        self.gmx.density('nvt.xtc', 'nvt.tpr', begin=when, silent=True)
        df = pd.read_table('density.xvg', skiprows=24, names=['Density'], index_col=0, sep='\s+')
        density_series = df.Density

        is_interface, bounds = is_interface(density_series)
        if not is_interface:
            raise Exception('No interface detected')
        z_liq_span = bounds[2] - bounds[1]
        z_liq = [bounds[1] + z_liq_span * 0.2, bounds[2] - z_liq_span * 0.2]
        z_gas_span_1 = bounds[1] - bounds[0]
        z_gas_span_2 = bounds[3] - bounds[2]
        z_gas_1 = [bounds[0], bounds[1] - z_gas_span_1 * 0.2]
        z_gas_2 = [bounds[2] + z_gas_span_2 * 0.2, bounds[3]]

        density_liq = density_series.loc[z_liq[0]:z_liq[1]].mean()
        density_gas = density_series.loc[z_gas_1[0]:z_gas_1[1]].append(density_series.loc[z_gas_2[0]:z_gas_2[1]]).mean()

        surface_tension, st_stderr = self.gmx.get_property_and_stderr('nvt.edr', '#Surf*SurfTen', begin=when)

        return {
            'simulation_length': potential_series.index[-1],
            'converged_from': when,
            'surface_tension': [surface_tension / 20, st_stderr / 20],
            'density_liq': density_liq / 1000,
            'density_gas': density_gas / 1000,
        }

    def clean(self):
        for f in os.listdir(os.getcwd()):
            if f.startswith('em.') or f.startswith('anneal.') or f.startswith('eq.') \
                    or f.endswith('.dfo') or f.startswith('#'):
                try:
                    os.remove(f)
                except:
                    pass
