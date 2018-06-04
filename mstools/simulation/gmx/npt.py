import os
import shutil

from .gmx import GmxSimulation
from ...analyzer import is_converged, block_average
from ...wrapper.ppf import delta_ppf


class Npt(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'npt'
        self.requirement = []
        self.logs = ['npt.log', 'hvap.log']
        self.n_atoms_default = 3000

    def build(self, export=True, ppf=None, minimize=False):
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, self.pdb, size=[i - 2 for i in self.box], silent=True)
        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd, mol_corr=self.pdb, size=self.box)

        # build msd for fast export
        self.packmol.build_box(self.pdb_list, [1] * len(self.pdb_list), self._single, size=self.box,
                               inp_file='build_single.inp', silent=True)
        self.dff.build_box_after_packmol(self.mol2_list, [1] * len(self.pdb_list), self._single,
                                         mol_corr=self._single, size=self.box)

        if export:
            self.fast_export_single(ppf=ppf, gro_out='_single.gro', top_out='topol.top')
            self.gmx.pdb2gro(self.pdb, 'conf.gro', [i / 10 for i in self.box], silent=True) # A to nm
            self.gmx.modify_top_mol_numbers('topol.top', self.n_mol_list)

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=298, P=1, jobname=None,
                dt=0.002, nst_eq=int(4E5), nst_run=int(5E5), nst_edr=100, nst_trr=int(5E4), nst_xtc=int(1E3),
                drde=False, **kwargs) -> [str]:
        if os.path.abspath(model_dir) != os.getcwd():
            shutil.copy(os.path.join(model_dir, gro), gro)
            shutil.copy(os.path.join(model_dir, top), top)
            for f in os.listdir(model_dir):
                if f.endswith('.itp'):
                    shutil.copy(os.path.join(model_dir, f), '.')

        if drde:
            ### Temperature dependent parameters
            # TODO Assumes ppf file named ff.ppf
            if os.path.abspath(model_dir) != os.getcwd():
                shutil.copy(os.path.join(model_dir, self._single), self._single)
            delta_ppf(os.path.join(model_dir, 'ff.ppf'), 'ff.ppf', T)
            mol_numbers = self.gmx.get_top_mol_numbers(top)
            self.fast_export_single(ppf='ff.ppf', gro_out='_single.gro', top_out=top)
            self.gmx.modify_top_mol_numbers(top, [n for m, n in mol_numbers])

        nprocs = self.jobmanager.nprocs
        commands = []
        # energy minimization
        self.gmx.prepare_mdp_from_template('t_em.mdp', mdp_out='grompp-em.mdp')
        cmd = self.gmx.grompp(mdp='grompp-em.mdp', gro=gro, top=top, tpr_out='em.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='em', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NVT annealing from 0 to 1000 K to target T with Langevin thermostat
        self.gmx.prepare_mdp_from_template('t_nvt_anneal.mdp', mdp_out='grompp-anneal.mdp', T=T,
                                           nsteps=int(1E5), nstxtcout=0)
        cmd = self.gmx.grompp(mdp='grompp-anneal.mdp', gro='em.gro', top=top, tpr_out='anneal.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='anneal', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NPT equilibrium with Langevin thermostat and Berendsen barostat
        self.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-eq.mdp', T=T, P=P,
                                           nsteps=nst_eq, nstxtcout=0, restart=True, pcoupl='berendsen')
        cmd = self.gmx.grompp(mdp='grompp-eq.mdp', gro='anneal.gro', top=top, tpr_out='eq.tpr',
                              cpt='anneal.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='eq', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NPT production with Langevin thermostat and Parrinello-Rahman barostat
        self.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-npt.mdp', T=T, P=P,
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
        # Use OpenMP instead of MPI when rerun hvap
        cmd = self.gmx.mdrun(name='hvap', nprocs=nprocs, n_omp=nprocs, rerun='npt.xtc', get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def extend(self, extend=500, jobname=None, sh=None) -> [str]:
        '''
        extend simulation for 500 ps
        '''
        self.gmx.extend_tpr('npt.tpr', extend, silent=True)

        nprocs = self.jobmanager.nprocs
        commands = []
        # Extending NPT production
        cmd = self.gmx.mdrun(name='npt', nprocs=nprocs, extend=True, get_cmd=True)
        commands.append(cmd)

        # Rerun enthalpy of vaporization
        commands.append('export GMX_MAXCONSTRWARN=-1')
        # Use OpenMP instead of MPI when rerun hvap
        cmd = self.gmx.mdrun(name='hvap', nprocs=nprocs, n_omp=nprocs, rerun='npt.xtc', get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure, sh=sh)
        return commands

    def analyze(self, check_converge=True, **kwargs):
        import numpy as np
        from ...panedr import edr_to_df

        df = edr_to_df('npt.edr')
        potential_series = df.Potential
        density_series = df.Density
        df = edr_to_df('hvap.edr')
        e_inter_series = (df.Potential)
        length = potential_series.index[-1]

        ### Check structure freezing using Diffusion. Only use last 400 ps data
        diffusion, _ = self.gmx.diffusion('npt.xtc', 'npt.tpr', begin=length - 400)
        if diffusion < 1E-8:  # cm^2/s
            return {
                'failed': True,
                'reason': 'freeze'
            }

        ### Check convergence
        if not check_converge:
            when = 0
        else:
            # use potential to do a initial determination
            # use at least 4/5 of the data
            _, when_pe = is_converged(potential_series, frac_min=0)
            when_pe = min(when_pe, length * 0.2)
            # use density to do a final determination
            _, when_dens = is_converged(density_series, frac_min=0)
            when = max(when_pe, when_dens)
            if when > length * 0.5:
                return None

        ### Get expansion and compressibility using fluctuation method
        nblock = 5
        blocksize = (length - when) / nblock
        expan_list = []
        compr_list = []
        for i in range(nblock):
            begin = when + blocksize * i
            end = when + blocksize * (i + 1)
            expan, compr = self.gmx.get_fluct_props('npt.edr', begin=begin, end=end)
            expan_list.append(expan)
            compr_list.append(compr)
        expansion, expan_stderr = np.mean(expan_list), np.std(expan_list, ddof=1) / np.sqrt(nblock)
        compressi, compr_stderr = np.mean(compr_list), np.std(compr_list, ddof=1) / np.sqrt(nblock)
        expan_stderr = float('%.1e' % expan_stderr)  # 2 effective number for stderr
        compr_stderr = float('%.1e' % compr_stderr)  # 2 effective number for stderr

        temperature_and_stderr, pressure_and_stderr, potential_and_stderr, density_and_stderr = \
            self.gmx.get_properties_stderr('npt.edr',
                                           ['Temperature', 'Pressure', 'Potential', 'Density'],
                                           begin=when)
        return {
            'simulation_length': length,
            'converged_from'   : when,
            'temperature'      : temperature_and_stderr,  # K
            'pressure'         : pressure_and_stderr,  # bar
            'potential'        : potential_and_stderr,  # kJ/mol
            'density'          : [i / 1000 for i in density_and_stderr],  # g/mL
            'e_inter'          : list(block_average(e_inter_series.loc[when:])),  # kJ/mol
            'expansion'        : [expansion, expan_stderr],
            'compressibility'  : [compressi, compr_stderr],
        }

    def clean(self):
        for f in os.listdir(os.getcwd()):
            if f.startswith('em.') or f.startswith('anneal.') or f.startswith('eq.') \
                    or f.endswith('.dfo') or f.startswith('#'):
                try:
                    os.remove(f)
                except:
                    pass

    @staticmethod
    def post_process(T_list, P_list, result_list, **kwargs) -> (dict, str):
        if len(set(T_list)) < 5 or len(set(P_list)) < 5:
            return None, 'T or P points less than 5'

        from mstools.analyzer.fitting import polyfit_2d

        density_list = []
        e_inter_list = []
        compres_list = []
        for i, result in enumerate(result_list):
            density_list.append(result['density'][0])
            e_inter_list.append(result['e_inter'][0])
            compres_list.append(result['compressibility'][0])

        coeff_density, score_density = polyfit_2d(T_list, P_list, density_list, 4)
        post_result = {'density-poly4-coeff': list(coeff_density),
                       'density-poly4-score': score_density,
                       }
        coeff_e_inter, score_e_inter = polyfit_2d(T_list, P_list, e_inter_list, 4)
        post_result.update({'e_inter-poly4-coeff': list(coeff_e_inter),
                            'e_inter-poly4-score': score_e_inter,
                            })

        return post_result, 'density-poly4-score %.4f e_inter-poly4-score %.4f' % (score_density, score_e_inter)

    @staticmethod
    def get_post_data(post_result, T, P, smiles_list, n_mol_list, **kwargs) -> dict:
        from mstools.analyzer.fitting import polyval_derivative_2d
        coeff_density = post_result['density-poly4-coeff']
        score_density = post_result['density-poly4-score']

        density, dDdT, dDdP = polyval_derivative_2d(T, P, 4, coeff_density)  # g/mL
        expansion = -1 / density * dDdT  # K^-1
        compressibility = 1 / density * dDdP  # bar^-1

        coeff_e_inter = post_result['e_inter-poly4-coeff']
        score_e_inter = post_result['e_inter-poly4-score']

        e_inter, dEdT, dEdP = polyval_derivative_2d(T, P, 4, coeff_e_inter)
        e_inter /= n_mol_list[0]  # kJ/mol
        hvap = 8.314 * T / 1000 - e_inter  # kJ/mol
        cp_inter = dEdT * 1000 / n_mol_list[0]  # J/mol.K

        import pybel
        py_mol = pybel.readstring('smi', smiles_list[0])
        cp_pv = - py_mol.molwt * P / density ** 2 * dDdT * 0.1  # J/mol/K

        return {
            'density'        : density,
            'expansion'      : expansion,
            'compressibility': compressibility,
            'e_inter'        : e_inter,
            'hvap'           : hvap,
            'cp_inter'       : cp_inter,
            'cp_pv'          : cp_pv,
            'score-density'  : score_density,
            'score-e_inter'  : score_e_inter,
        }
