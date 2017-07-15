import os
import shutil

from mstools.analyzer.series import is_converged
from .gmx import GmxSimulation


class NvtVacuum(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'nvt-vacuum'
        self.logs = ['nvt-vacuum.log']
        self.requirement = []

    def build(self, minimize=False):
        # to calculate enthalpy of vaporization, only one molecule in vacuum of length 30 A
        self.n_mol_list = [1 for i in self.n_mol_list]
        self.length = 30
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, 'init.pdb', length=self.length - 2, silent=True)

        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd, mol_corr='init.pdb',
                                         length=self.length)
        self.export(minimize=True, vacuum=False)

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=None, P=None, jobname=None, **kwargs):
        if os.path.abspath(model_dir) != os.getcwd():
            shutil.copy(os.path.join(model_dir, gro), gro)
            shutil.copy(os.path.join(model_dir, top), top)
            for f in os.listdir(model_dir):
                if f.endswith('.itp'):
                    shutil.copy(os.path.join(model_dir, f), '.')

        commands = []
        self.gmx.prepare_mdp_from_template('t_nvt_vacuum.mdp', T=T, nsteps=int(1E6))
        cmd = self.gmx.grompp(gro=gro, top=top, tpr_out=self.procedure, get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name=self.procedure, nprocs=1, get_cmd=True)
        commands.append(cmd)
        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)

    def analyze(self, dirs=None):
        if dirs is None:
            dirs = ['.']
        import panedr
        import pandas as pd
        import numpy as np

        # TODO check convergence, utilizing previous cycles
        temp_series = pd.Series()
        pe_series = pd.Series()
        for dir in dirs:
            df = panedr.edr_to_df(os.path.join(dir, self.procedure + '.edr'))
            temp_series = temp_series.append(df.Temperature)
            pe_series = pe_series.append(df.Potential)

        converged, when = is_converged(pe_series)

        return converged, {
            'temperature': np.mean(temp_series[when:]),
            'potential': np.mean(pe_series[when:]),
        }
