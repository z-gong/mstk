import os
import shutil

from .gmx import GmxSimulation
from ...analyzer import is_converged, block_average
from ...wrapper.ppf import delta_ppf


class Npt(GmxSimulation):
    '''
    Npt simulation with GROMACS enables the calculation of equation of states, cohesive energy,
    which can be used to derive heat of vaporization, thermal expansion, compressibility, etc...

    Parameters
    ----------
    packmol : Packmol
    dff : DFF
    gmx : GMX
    jobmanager : subclass of JobManager

    Attributes
    ----------
    packmol : Packmol
    dff : Packmol
    gmx : GMX
    jobmanager : subclass of JobManager
    n_atom_default : int
    n_mol_default : int
    n_mol_list : list of int
    logs : list of str
    '''

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'npt'
        self.logs = ['npt.log', 'hvap.log']
        self.n_atom_default = 3000
        self.n_mol_default = 75

    def set_system(self, smiles_list, **kwargs):
        super().set_system(smiles_list, **kwargs)

    def build(self, export=True, ppf=None):
        '''
        Generate control scripts for simulation engine. e.g. mdp files for GROMACS.
        Also generate force field files for target temperatures if force field is temperature dependent.

        Parameters
        ----------
        export
        ppf

        Returns
        -------

        '''
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, self.pdb,
                               size=[i - 2 for i in self.box], silent=True)
        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd,
                                         mol_corr=self.pdb, size=self.box)

        # build msd for fast export
        self.packmol.build_box(self.pdb_list, [1] * len(self.pdb_list), self._single_pdb,
                               size=self.box,
                               inp_file='build_single.inp', silent=True)
        self.dff.build_box_after_packmol(self.mol2_list, [1] * len(self.pdb_list), self._single_msd,
                                         mol_corr=self._single_pdb, size=self.box)

        if export:
            self.fast_export_single(ppf=ppf, gro_out='_single.gro', top_out='topol.top')
            self.gmx.pdb2gro(self.pdb, 'conf.gro', [i / 10 for i in self.box],
                             silent=True)  # A to nm
            self.gmx.modify_top_mol_numbers('topol.top', self.n_mol_list)

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=298, P=1, jobname=None,
                TANNEAL=800,
                dt=0.002, nst_eq=int(4E5), nst_run=int(5E5), nst_edr=100, nst_trr=int(5E4),
                nst_xtc=int(1E3),
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
                shutil.copy(os.path.join(model_dir, self._single_msd), self._single_msd)
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

        gro_em = 'em.gro'
        # NVT annealing from 0 to TANNEAL to target T with Langevin thermostat
        if TANNEAL is not None:
            self.gmx.prepare_mdp_from_template('t_nvt_anneal.mdp', mdp_out='grompp-anneal.mdp', T=T,
                                               TANNEAL=TANNEAL,
                                               nsteps=int(1E5), nstxtcout=0)
            cmd = self.gmx.grompp(mdp='grompp-anneal.mdp', gro='em.gro', top=top,
                                  tpr_out='anneal.tpr', get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name='anneal', nprocs=nprocs, get_cmd=True)
            commands.append(cmd)

            gro_em = 'anneal.gro'

        # NPT equilibrium with Langevin thermostat and Berendsen barostat
        self.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-eq.mdp', T=T, P=P,
                                           nsteps=nst_eq, nstxtcout=0, pcoupl='berendsen')
        cmd = self.gmx.grompp(mdp='grompp-eq.mdp', gro=gro_em, top=top, tpr_out='eq.tpr',
                              get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='eq', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NPT production with Langevin thermostat and Parrinello-Rahman barostat
        self.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-npt.mdp', T=T, P=P,
                                           dt=dt, nsteps=nst_run, nstenergy=nst_edr,
                                           nstxout=nst_trr, nstvout=nst_trr,
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
        self.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-hvap.mdp', nstxtcout=0,
                                           restart=True)
        cmd = self.gmx.grompp(mdp='grompp-hvap.mdp', gro='eq.gro', top=top_hvap, tpr_out='hvap.tpr',
                              get_cmd=True)
        commands.append(cmd)
        # Use OpenMP instead of MPI when rerun hvap
        cmd = self.gmx.mdrun(name='hvap', nprocs=nprocs, n_omp=nprocs, rerun='npt.xtc',
                             get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def run(self):
        super().run()

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
        cmd = self.gmx.mdrun(name='hvap', nprocs=nprocs, n_omp=nprocs, rerun='npt.xtc',
                             get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure, sh=sh)
        return commands

    def analyze(self, check_converge=True, **kwargs):
        import numpy as np
        from ...panedr import edr_to_df

        df = edr_to_df('npt.edr')
        potential_series = df.Potential
        density_series = df.Density
        t_real = df.Temperature.mean()
        length = potential_series.index[-1]

        df_hvap = edr_to_df('hvap.edr')
        einter_series = (df_hvap.Potential)

        ### Check the ensemble. KS test on the distribution of kinetic energy
        ### TODO This can be optimized
        import physical_validation as pv
        parser = pv.data.GromacsParser(self.gmx.GMX_BIN)
        data = parser.get_simulation_data(mdp='grompp-npt.mdp', top='topol.top', edr='npt.edr')
        p = pv.kinetic_energy.distribution(data, strict=True, verbosity=0)
        # If test does not pass, set the desired temperature to t_real.
        # Because small deviation in temperature exists for Langevin thermostat
        # 3 Kelvin of deviation is permitted
        if p < 0.01 and abs(data.ensemble.temperature - t_real) < 3:
            try:
                data._SimulationData__ensemble._EnsembleData__t = t_real
                p = pv.kinetic_energy.distribution(data, strict=True, verbosity=0)
            except Exception as e:
                print(repr(e))
        if p < 0.01:
            return {
                'failed': True,
                'reason': 'incorrect_ensemble'
            }

        ### Check structure freezing using Density
        if density_series.min() / 1000 < 0.1:  # g/mL
            return {
                'failed': True,
                'reason': 'vaporize'
            }

        ### Check structure freezing using Diffusion of COM of molecules. Only use last 400 ps data
        diffusion, _ = self.gmx.diffusion('npt.xtc', 'npt.tpr', mol=True, begin=length - 400)
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
            'length'     : length,
            'converge'   : when,
            'temperature': temperature_and_stderr,  # K
            'pressure'   : pressure_and_stderr,  # bar
            'potential'  : potential_and_stderr,  # kJ/mol
            'density'    : [i / 1000 for i in density_and_stderr],  # g/mL
            'einter'     : list(block_average(einter_series.loc[when:])),  # kJ/mol
            'expansion'  : [expansion, expan_stderr],
            'compress'   : [compressi, compr_stderr],
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
    def post_process(T_list, P_list, result_list, n_mol_list, **kwargs) -> (dict, str):
        t_set = set(T_list)
        p_set = set(P_list)
        if len(t_set) < 5 or len(p_set) < 5:
            return None, 'T or P points less than 5'

        from mstools.analyzer.fitting import polyfit_2d, polyfit

        def round3(x):
            return float('%.3e' % x)

        ### einter divided by number of molecules
        dens_stderr_list = [list(map(round3, result['density'])) for result in result_list]
        eint_stderr_list = [list(map(lambda x: round3(x / n_mol_list[0]), result['einter'])) for
                            result in result_list]
        comp_stderr_list = [list(map(round3, result['compress'])) for result in result_list]

        dens_list = [i[0] for i in dens_stderr_list]
        eint_list = [i[0] for i in eint_stderr_list]
        comp_list = [i[0] for i in comp_stderr_list]

        ### Fit against T, P with poly4
        coeff_dens, score_dens = polyfit_2d(T_list, P_list, dens_list, 4)
        coeff_eint, score_eint = polyfit_2d(T_list, P_list, eint_list, 4)

        ### Fit against only T with poly3
        t_p_dens_list = list(map(list, zip(T_list, P_list, dens_list)))
        t_p_dens_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
        t_p_eint_list = list(map(list, zip(T_list, P_list, eint_list)))
        t_p_eint_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
        t_p_comp_list = list(map(list, zip(T_list, P_list, comp_list)))
        t_p_comp_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T

        t_p_dens_stderr_list = list(map(list, zip(T_list, P_list, dens_stderr_list)))
        t_p_dens_stderr_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
        t_p_eint_stderr_list = list(map(list, zip(T_list, P_list, eint_stderr_list)))
        t_p_eint_stderr_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
        t_p_comp_stderr_list = list(map(list, zip(T_list, P_list, comp_stderr_list)))
        t_p_comp_stderr_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T

        t_dens_poly3 = {}
        t_eint_poly3 = {}
        t_comp_poly3 = {}

        for p in sorted(p_set):
            _t_list = [element[0] for element in t_p_dens_list if element[1] == p]
            _dens_list = [element[2] for element in t_p_dens_list if element[1] == p]
            _eint_list = [element[2] for element in t_p_eint_list if element[1] == p]
            _comp_list = [element[2] for element in t_p_comp_list if element[1] == p]

            if len(_t_list) < 5:
                continue

            # density-T relation is fitted by 3-order polynomial function
            _t_dens_coeff, _t_dens_score = polyfit(_t_list, _dens_list, 3)
            _t_eint_coeff, _t_eint_score = polyfit(_t_list, _eint_list, 3)
            _t_comp_coeff, _t_comp_score = polyfit(_t_list, _comp_list, 3)

            t_dens_poly3[p] = [list(map(round3, _t_dens_coeff)), round3(_t_dens_score),
                               min(_t_list), max(_t_list)]
            t_eint_poly3[p] = [list(map(round3, _t_eint_coeff)), round3(_t_eint_score),
                               min(_t_list), max(_t_list)]
            t_comp_poly3[p] = [list(map(round3, _t_comp_coeff)), round3(_t_comp_score),
                               min(_t_list), max(_t_list)]

        post_result = {
            'density'         : t_p_dens_stderr_list,
            'einter'          : t_p_eint_stderr_list,
            'compress'        : t_p_comp_stderr_list,
            'density-poly4'   : [list(map(round3, coeff_dens)), round3(score_dens)],
            'einter-poly4'    : [list(map(round3, coeff_eint)), round3(score_eint)],
            'density-t-poly3' : t_dens_poly3,
            'einter-t-poly3'  : t_eint_poly3,
            'compress-t-poly3': t_comp_poly3,
        }

        return post_result, 'density-poly4-score %.4f einter-poly4-score %.4f' % (
            score_dens, score_eint)

    @staticmethod
    def get_post_data(post_result, T, P, smiles_list, **kwargs) -> dict:
        from mstools.analyzer.fitting import polyval_derivative_2d, polyval, polyval_derivative, \
            polyfit

        ### Calculate with T,P-poly4. Not accurate enough, especially for expansion and compressibility
        coeff_dens, score_dens = post_result['density-poly4']
        coeff_eint, score_eint = post_result['einter-poly4']
        density4, dDdT4, dDdP4 = polyval_derivative_2d(T, P, 4, coeff_dens)  # g/mL
        einter4, dEdT4, dEdP4 = polyval_derivative_2d(T, P, 4, coeff_eint)  # kJ/mol
        expansion4 = -1 / density4 * dDdT4  # K^-1
        compressibility4 = 1 / density4 * dDdP4  # bar^-1

        import pybel
        py_mol = pybel.readstring('smi', smiles_list[0])
        cp_inter4 = dEdT4 * 1000  # J/mol.K
        cp_pv4 = - py_mol.molwt * P / density4 ** 2 * dDdT4 * 0.1  # J/mol/K

        ### T-poly3
        _p_dens_list = []
        _p_eint_list = []
        _p_comp_list = []
        _p_dDdT_list = []
        _p_dEdT_list = []
        for _p in post_result['density-t-poly3']:
            coef, score, tmin, tmax = post_result['density-t-poly3'][str(_p)]
            if score < 0.999 or T < tmin - 10 or T > tmax + 10:
                continue

            dens, dDdT = polyval_derivative(T, coef)
            _p_dens_list.append([float(_p), dens])
            _p_dDdT_list.append([float(_p), dDdT])

        for _p in post_result['einter-t-poly3']:
            coef, score, tmin, tmax = post_result['einter-t-poly3'][str(_p)]
            if score < 0.999 or T < tmin - 10 or T > tmax + 10:
                continue

            eint, dEdT = polyval_derivative(T, coef)
            _p_eint_list.append([float(_p), eint])
            _p_dEdT_list.append([float(_p), dEdT])

        for _p in post_result['compress-t-poly3']:
            coef, score, tmin, tmax = post_result['compress-t-poly3'][str(_p)]
            if score < 0.95 or T < tmin - 10 or T > tmax + 10:
                continue

            _p_comp_list.append([float(_p), polyval(T, coef)])

        ### Default value
        density = None
        einter = None
        hvap = None
        cp_inter = None
        cp_pv = None
        expansion = None
        compressibility = None

        if len(_p_dens_list) >= 5:
            coef, score = polyfit(*zip(*_p_dens_list), 3)
            _p_list = list(zip(*_p_dens_list))[0]
            if P > min(_p_list) - 10 and P < max(_p_list) + 10:
                density = polyval(P, coef)

                coef, score = polyfit(*zip(*_p_dDdT_list), 3)
                dDdT = polyval(P, coef)
                expansion = -1 / density * dDdT  # K^-1
                cp_pv = - py_mol.molwt * P / density ** 2 * dDdT * 0.1  # J/mol/K

        if len(_p_eint_list) >= 5:
            coef, score = polyfit(*zip(*_p_eint_list), 3)
            _p_list = list(zip(*_p_eint_list))[0]
            if P > min(_p_list) - 10 and P < max(_p_list) + 10:
                einter = polyval(P, coef)
                hvap = 8.314 * T / 1000 - einter  # kJ/mol

                coef, score = polyfit(*zip(*_p_dEdT_list), 3)
                dEdT = polyval(P, coef)
                cp_inter = dEdT * 1000  # J/mol.K

        if len(_p_comp_list) >= 5:
            coef, score = polyfit(*zip(*_p_comp_list), 3)
            _p_list = list(zip(*_p_comp_list))[0]
            if P > min(_p_list) - 10 and P < max(_p_list) + 10:
                compressibility = polyval(P, coef)

        return {
            'density'            : density,
            'einter'             : einter,
            'hvap'               : hvap,
            'cp_inter'           : cp_inter,
            'cp_pv'              : cp_pv,
            'expansion'          : expansion,
            'compress'           : compressibility,
            'density-poly4-score': score_dens,
            'einter-poly4-score' : score_eint,
            'density-poly4'      : density4,
            'einter-poly4'       : einter4,
            'cp_inter-poly4'     : cp_inter4,
            'cp_pv-poly4'        : cp_pv4,
            'expansion-poly4'    : expansion4,
            'compress-poly4'     : compressibility4,
        }
