import os
import shutil

from .gmx import GmxSimulation
from ...analyzer.series import block_average
from ...analyzer.structure import N_vaporize_condense
from ...analyzer.fitting import *
from ...wrapper.ppf import delta_ppf


class NvtSlab(GmxSimulation):
    '''
    Nvt slab simulation with GROMACS enables the calculation of surface tension and liquid-vapor phase equilibrium,
    which can be used to derive critical temperature and critical density.

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
        self.procedure = 'nvt-slab'
        self.logs = ['nvt.log']
        self.n_atom_default = 12000
        self.n_mol_default = 400

    def set_system(self, smiles_list, **kwargs):
        super().set_system(smiles_list, **kwargs)

    def build(self, export=True, ppf=None, length=None):
        if length is None:
            ### Slim box (z=2.0*x) in case of fully vaporization at high T
            length = (self.vol / 2.0) ** (1 / 3)
        self.box = [length, length, length * 6.0]
        slab = self.vol / length ** 2

        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, 'init.pdb',
                               size=[i - 2 for i in self.box], slab=slab, silent=True)

        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd,
                                         mol_corr='init.pdb', size=self.box)

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

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=298, jobname=None,
                TANNEAL=None,
                dt=0.002, nst_eq=int(4E5), nst_run=int(4E6), nst_edr=100, nst_trr=int(5E4),
                nst_xtc=int(5E2),
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

        # NVT equilibrium
        self.gmx.prepare_mdp_from_template('t_nvt_slab.mdp', mdp_out='grompp-eq.mdp', T=T,
                                           nsteps=nst_eq, nstxtcout=0)
        cmd = self.gmx.grompp(mdp='grompp-eq.mdp', gro=gro_em, top=top, tpr_out='eq.tpr',
                              get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='eq', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NVT production
        self.gmx.prepare_mdp_from_template('t_nvt_slab.mdp', mdp_out='grompp-nvt.mdp', T=T,
                                           dt=dt, nsteps=nst_run, nstenergy=nst_edr,
                                           nstxout=nst_trr, nstvout=nst_trr,
                                           nstxtcout=nst_xtc, restart=True)
        cmd = self.gmx.grompp(mdp='grompp-nvt.mdp', gro='eq.gro', top=top, tpr_out='nvt.tpr',
                              cpt='eq.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='nvt', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def run(self):
        super().run()

    def extend(self, extend=2000, jobname=None, sh=None) -> [str]:
        '''
        extend simulation for 2000 ps
        '''
        self.gmx.extend_tpr('nvt.tpr', extend, silent=True)

        nprocs = self.jobmanager.nprocs
        commands = []
        # Extending NVT production
        cmd = self.gmx.mdrun(name='nvt', nprocs=nprocs, extend=True, get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure, sh=sh)
        return commands

    def analyze(self, debug=False, **kwargs):
        import pandas as pd
        from ...panedr import edr_to_df
        from ...analyzer.series import is_converged
        from ...analyzer.structure import check_vle_density

        df = edr_to_df('nvt.edr')
        potential_series = df.Potential
        length = potential_series.index[-1]

        ### Check structure freezing using Diffusion of COM of molecules. Only use last 400 ps of data
        diffusion, _ = self.gmx.diffusion('nvt.xtc', 'nvt.tpr', mol=True, begin=length - 400)
        if diffusion < 1E-8:  # cm^2/s
            return {
                'failed': True,
                'reason': 'freeze'
            }

        # use potential to do a initial determination
        # use at least 4/5 of the data
        _, when = is_converged(potential_series, frac_min=0)
        when = min(when, length * 0.2)

        # cut trj to pieces, analyze the change of density of liq and gas phases
        dt = potential_series.index[1] - potential_series.index[0]
        n_frame = len(potential_series.loc[when:])
        dt_piece = 200  # ps
        # use int in stead of math.ceil. make sure each peace has at least dt_piece
        n_pieces = int((potential_series.index[-1] - when) / dt_piece)

        if debug:
            print('n_frame=%i dt=%f dt_piece=%f n_pieces=%i' % (n_frame, dt, dt_piece, n_pieces))
            self.lz_liq_series = pd.Series()
            self.lz_gas_series = pd.Series()
            self.lz_int_series = pd.Series()

            import mdtraj
            time: [float] = []  # time of frames, saved to calculate vaporize
            Ngas: [int] = []  # number of molecules in gas phases
            phases: [[str]] = []  # [['l', 'g','i'], ...] # n_atoms * n_frames

            def in_which_phase(z):
                z -= idxmax
                while z < 0:
                    z += lz
                while z > lz:
                    z -= lz
                if z < r1 - s1 * mul or z > r2 + s2 * mul:
                    return 'l'  # liquid phase
                elif z > r1 + s1 * mul + 1 and z < r2 - s2 * mul - 1:  # shrink 1 nm for gas phase near interface
                    return 'g'  # gas phase
                else:
                    return 'i'  # interface

            def in_gas(z):
                z -= idxmax
                while z < 0:
                    z += lz
                while z > lz:
                    z -= lz
                if z > r1 + s1 * mul and z < r2 - s2 * mul:
                    return True
                return False

        dliq_series = pd.Series()  # density of liquid phase timeseries
        dgas_series = pd.Series()  # density of liquid phase timeseries
        for n in range(n_pieces):
            xvg = 'density-%i.xvg' % n
            begin = when + dt_piece * n
            end = when + dt_piece * (n + 1) - dt
            if debug:
                print('Time: ', begin, end)

            self.gmx.density('nvt.xtc', 'nvt.tpr', xvg=xvg, begin=begin, end=end, silent=True)
            df = pd.read_table(xvg, skiprows=24, names=['Density'], index_col=0, sep='\s+')
            density_series = df.Density
            if not debug:
                os.remove(xvg)

            dz = density_series.index[1] - density_series.index[0]
            lz = dz * len(density_series)

            # Move gas phase to center. Series index starts from 0
            idxmax = density_series.idxmax()
            density_series.index -= density_series.idxmax()
            _index_list = density_series.index.tolist()
            for i, v in enumerate(_index_list):
                if v < 0:
                    _index_list[i] += lz
            density_series.index = _index_list
            density_series = density_series.sort_index()

            is_interface, is_gas_center, nodes = check_vle_density(density_series)
            if not is_interface or not is_gas_center:
                continue

            # Move gas phase to exact center. Series index not starts from 0
            center = sum(nodes) / 2
            _index_list = density_series.index.tolist()
            if center < lz / 2:
                shift = lz / 2 - center
                for i, v in enumerate(_index_list):
                    if v > lz - shift:
                        _index_list[i] -= lz
            else:
                shift = center - lz / 2
                for i, v in enumerate(_index_list):
                    if v < shift:
                        _index_list[i] += lz
            density_series.index = _index_list
            density_series = density_series.sort_index()

            # Fit Tanh to determine the thickness of liquid phase
            dens_series_left = density_series.loc[:center]
            dens_series_righ = density_series.loc[center:]
            _max = max(density_series)
            _min = min(density_series)
            _c = (_max + _min) / 2
            _A = (_max - _min) / 2
            _coef1, _score1 = fit_vle_tanh(dens_series_left.index, dens_series_left,
                                           guess=[_c, -_A, nodes[0], 1])
            _coef2, _score1 = fit_vle_tanh(dens_series_righ.index, dens_series_righ,
                                           guess=[_c, _A, nodes[1], 1])
            # if debug:
            #     print(dens_series_left)
            #     print(dens_series_righ)
            #     print('Tanh: ', _coef1)
            #     print('Tanh: ', _coef2)
            c1, A1, r1, s1 = _coef1  # A1 is negative
            c2, A2, r2, s2 = _coef2  # A2 is positive

            # Check if density fluctuate along the z axis
            if abs(c1 - A1 - density_series.max()) > 50 / 1000 * density_series.max() or \
                    abs(c2 + A2 - density_series.max()) > 50 / 1000 * density_series.max() or \
                    abs(c1 + A1 - density_series.min()) > 100 / 1000 * density_series.max() or \
                    abs(c2 - A2 - density_series.min()) > 100 / 1000 * density_series.max():
                continue

            mul = 2.6466524123622457  # arctanh(0.99)
            lz_liq = r1 - density_series.index[0] + density_series.index[
                -1] - r2 + dz - s1 * mul - s2 * mul
            lz_gas = r2 - r1 - s1 * mul - s2 * mul
            lz_int = (s1 + s2) * mul * 2  # thickness of two interfaces
            if debug:
                print('Thickness of liquid, gas and interfaces: ', lz_liq, lz_gas, lz_int)
                self.lz_liq_series.at[begin] = lz_liq
                self.lz_gas_series.at[begin] = lz_gas
                self.lz_int_series.at[begin] = lz_int

            # Liquid phase should be at least 2 nm. Gas phase should be at least 5 nm
            if lz_liq < 2 or lz_gas < 5:
                continue

            d_liq = (c1 - A1 + c2 + A2) / 2
            d_gas = (c1 + A1 + c2 - A2) / 2
            d_gas = max(d_gas, 0)

            dliq_series.at[begin] = d_liq
            dgas_series.at[begin] = d_gas

            if debug:
                gro_com = 'com-%i.gro' % n
                self.gmx.traj_com('nvt.xtc', 'nvt.tpr', gro_com, begin=begin, end=end, silent=True)
                trj: mdtraj.Trajectory = mdtraj.load(gro_com)

                if phases == []:
                    phases = [[] for i in range(trj.n_atoms)]
                for n_frame, xyz_array in enumerate(trj.xyz):
                    time.append(trj.time[n_frame])
                    Ngas.append(0)
                    for n_atom, xyz in enumerate(xyz_array):
                        phase = in_which_phase(xyz[2])
                        phases[n_atom].append(phase)
                        if in_gas(xyz[2]):
                            Ngas[-1] += 1
                        # if n_frame == 0:
                        #     print(xyz, phase, idxmax, lz, nodes, lz_interface)
                    # import sys
                    # sys.exit()

        if debug:
            self.dliq_series = dliq_series
            self.dgas_series = dgas_series
            self.Ngas_series = pd.Series(Ngas, time)
            self.phase_series_list = [pd.Series(l, time) for l in phases]
            print(dliq_series)
            print(dgas_series)

        # Failed if more than 1/4 pieces do not have interface
        # Failed if pieces at last 1/5 time span have no interface
        if len(dliq_series) < n_pieces * 0.75 or dliq_series.index[-1] < length * 0.8:
            return {
                'failed': True,
                'reason': 'no_interface'
            }

        _, when_liq = is_converged(dliq_series, frac_min=0)
        _, when_gas = is_converged(dgas_series, frac_min=0)
        # Convergence should be at least 4 ns
        if when_liq > length - 4000:
            return None
        if when_gas > length - 4000:
            if dgas_series.loc[when_gas:].mean() > 5:
                return None
            else:
                # Even if not converge. The density is so small < 5 kg/m^3. Considered as converged.
                when_gas = length - 4000
        when = max(when_liq, when_gas)

        if debug:
            self.N_vaporize = self.N_condense = 0
            for phase_series in self.phase_series_list:
                n_vap, n_con = N_vaporize_condense(phase_series.loc[when:])
                self.N_vaporize += n_vap
                self.N_condense += n_con

        temperature_and_stderr, pressure_and_stderr, pzz_and_stderr, st_and_stderr, potential_and_stderr = \
            self.gmx.get_properties_stderr('nvt.edr',
                                           ['Temperature', 'Pressure', 'Pres-ZZ', '#Surf*SurfTen',
                                            'Potential'],
                                           begin=when)
        return {
            'length'     : length,
            'converge'   : when,
            'temperature': temperature_and_stderr,  # K
            'pressure'   : pressure_and_stderr,  # bar
            'pzz'        : pzz_and_stderr,  # bar
            'st'         : [i / 20 for i in st_and_stderr],  # mN/m
            'potential'  : potential_and_stderr,  # kJ/mol
            'dliq'       : list(block_average(dliq_series.loc[when:] / 1000)),  # g/mL
            'dgas'       : list(block_average(dgas_series.loc[when:] / 1000)),  # g/mL
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
    def post_process(T_list, result_list, **kwargs) -> (dict, str):
        if len(set(T_list)) < 5:
            return None, 'T points less than 5'

        import numpy as np

        def round3(x):
            return float('%.3e' % x)

        dliq_stderr_list = [list(map(round3, result['dliq'])) for result in result_list]
        dgas_stderr_list = [list(map(round3, result['dgas'])) for result in result_list]
        st___stderr_list = [list(map(round3, result['st'])) for result in result_list]
        pzz__stderr_list = [list(map(round3, result['pzz'])) for result in result_list]

        dliq_list = [i[0] for i in dliq_stderr_list]
        dgas_list = [i[0] for i in dgas_stderr_list]
        st___list = [i[0] for i in st___stderr_list]
        pzz__list = [i[0] for i in pzz__stderr_list]

        coeff_dminus, score_dminus = fit_vle_dminus(T_list, np.array(dliq_list) - np.array(
            dgas_list))  # Tc, B
        Tc = coeff_dminus[0]
        if Tc < 0 or Tc > 1200:
            return None, 'Fit VLE density failed'

        coeff_dplus, score_dplus = fit_vle_dplus(T_list, np.array(dliq_list) + np.array(dgas_list),
                                                 Tc)  # Dc, A
        Dc = coeff_dplus[0]
        if Dc < 0 or Dc > 3.0:
            return None, 'Fit VLE density failed'

        coeff_st, score_st = fit_vle_st(T_list, st___list, Tc)  # Ast, n

        post_result = {
            'dliq'      : list(map(list, zip(T_list, dliq_stderr_list))),
            'dgas'      : list(map(list, zip(T_list, dgas_stderr_list))),
            'st'        : list(map(list, zip(T_list, st___stderr_list))),
            'pzz'       : list(map(list, zip(T_list, pzz__stderr_list))),
            'dminus-fit': [list(map(round3, coeff_dminus)), round3(score_dminus)],
            'dplus-fit' : [list(map(round3, coeff_dplus)), round3(score_dplus)],
            'st-fit'    : [list(map(round3, coeff_st)), round3(score_st)]
        }

        return post_result, 'Tc %.1f Dc %.3f' % (Tc, Dc)

    @staticmethod
    def get_post_data(post_result, T, **kwargs) -> dict:
        coeff_dminus, score_dminus = post_result['dminus-fit']
        coeff_dplus, score_dplus = post_result['dplus-fit']
        coeff_st, score_st = post_result['st-fit']

        Tc, _B = coeff_dminus
        Dc, _A = coeff_dplus

        if T > Tc:
            return {
                'tc'   : Tc,  # K
                'dc'   : Dc,  # g/mL
                'error': 'T larger than Tc',
            }

        dminus = vle_dminus(T, Tc, _B)
        dplus = vle_dplus(T, Dc, _A, Tc)
        dliq = (dplus + dminus) / 2
        dgas = (dplus - dminus) / 2

        _Ast, _n = coeff_st
        st = vle_st(T, _Ast, _n, Tc)

        return {
            'tc'  : Tc,  # K
            'dc'  : Dc,  # g/mL
            'dliq': dliq,  # g/mL
            'dgas': dgas,  # g/mL
            'st'  : st,  # mN/m
        }
