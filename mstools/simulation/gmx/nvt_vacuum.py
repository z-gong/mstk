import os
import shutil

from mstools.analyzer.series import is_converged
from .gmx import GmxSimulation
from ...wrapper.ppf import delta_ppf


class NvtVacuum(GmxSimulation):
    '''
    Nvt vacuum simulation with GROMACS enables the calculation of intramolecular enthalpy of single molecule.

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
        self.procedure = 'nvt-vacuum'
        self.logs = ['vacuum.log']
        self.n_atoms_default = 1

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

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=298, jobname=None,
                dt=0.002, nst_eq=int(4E5), nst_run=int(1E6), nst_edr=100, nst_trr=100, nst_xtc=0,
                drde=False, **kwargs):
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

        commands.append('export GMX_MAXCONSTRWARN=-1')
        # Energy minimization
        self.gmx.prepare_mdp_from_template('t_em_vacuum.mdp', mdp_out='grompp-em.mdp')
        cmd = self.gmx.grompp(mdp='grompp-em.mdp', gro=gro, top=top, tpr_out='em.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='em', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NVT eq
        self.gmx.prepare_mdp_from_template('t_nvt_vacuum.mdp', mdp_out='grompp-eq.mdp', T=T,
                                           nsteps=nst_eq, nstxtcout=0)
        cmd = self.gmx.grompp(mdp='grompp-eq.mdp', gro='em.gro', top=top, tpr_out='eq.tpr',
                              get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='eq', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NVT production
        self.gmx.prepare_mdp_from_template('t_nvt_vacuum.mdp', mdp_out='grompp-vacuum.mdp', T=T,
                                           dt=dt, nsteps=nst_run, nstenergy=nst_edr,
                                           nstxout=nst_trr, nstvout=nst_trr,
                                           nstxtcout=nst_xtc, restart=True)
        cmd = self.gmx.grompp(mdp='grompp-vacuum.mdp', gro='eq.gro', top=top, tpr_out='vacuum.tpr',
                              cpt='eq.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='vacuum', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def run(self):
        super().run()

    def analyze(self, dirs=None):
        if dirs is None:
            dirs = ['.']
        from ...panedr import edr_to_df
        import pandas as pd
        import numpy as np

        # TODO check convergence, utilizing previous cycles
        temp_series = pd.Series()
        pe_series = pd.Series()
        for dir in dirs:
            df = edr_to_df(os.path.join(dir, self.procedure + '.edr'))
            temp_series = temp_series.append(df.Temperature)
            pe_series = pe_series.append(df.Potential)

        converged, when = is_converged(pe_series)

        return converged, {
            'temperature': np.mean(temp_series[when:]),
            'potential'  : np.mean(pe_series[when:]),
        }
