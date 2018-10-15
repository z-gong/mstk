import os
import math

from .gauss import GaussSimulation
from ...utils import create_mol_from_smiles, generate_conformers


class Cv(GaussSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'cv'
        self.requirement = []
        self.logs = ['conf-0.log']

    def build(self, export=True, ppf=None, minimize=False):
        pass

    def prepare(self, gjf_name=None, n_conformer=0, T_list: [float] = None, jobname=None) -> [str]:
        if gjf_name is None:
            gjf_name = 'conf'
        self.logs = ['%s-%i.log' % (gjf_name, i) for i in range(n_conformer)]
        commands = []
        nprocs = self.jobmanager.nprocs

        # TODO Only first molecule
        smiles = self.smiles_list[0]
        py_mol = create_mol_from_smiles(smiles)

        if n_conformer == 0:
            name = '%s-%i' % (gjf_name, 0)
            self.gauss.generate_gjf_cv(name, py_mol, T_list=T_list)
            cmds = self.gauss.run_gjf(name + '.gjf', nprocs=nprocs, get_cmd=True)
            commands += cmds
        else:
            conformers = generate_conformers(py_mol, n_conformer)
            for i, conformer in enumerate(conformers):
                name = '%s-%i' % (gjf_name, i)
                self.gauss.generate_gjf_cv(name, conformer, T_list=T_list)
                cmds = self.gauss.run_gjf(name + '.gjf', nprocs=nprocs, get_cmd=True)
                commands += cmds

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def analyze(self, dirs=None, logs: [str] = None):
        def process_log(log):
            T_list = []
            scale_list = []
            Cv_list = []
            Cv_corr_list = []

            content = open(log).read()
            if content.find('Normal termination') == -1:
                return None
            if content.find('Error termination') > -1:
                return None
            if content.find('imaginary frequencies') > -1:
                return None

            f = open(log)
            while True:
                line = f.readline()
                if line == '':
                    break

                if line.strip().startswith('- Thermochemistry -'):
                    line = f.readline()
                    line = f.readline()
                    T = float(line.strip().split()[1])
                    T_list.append(T)
                    line = f.readline()
                    if line.strip().startswith('Thermochemistry will use frequencies scaled by'):
                        scale = float(line.strip().split()[-1][:-1])
                    else:
                        scale = 1
                    scale_list.append(scale)

                elif line.strip().startswith('E (Thermal)             CV                S'):
                    line = f.readline()
                    line = f.readline()
                    Cv = float(line.strip().split()[2]) * 4.184
                    Cv_list.append(Cv)
                    line = f.readline()
                    if line.strip().startswith('Corrected for'):
                        line = f.readline()
                        Cv_corr = float(line.strip().split()[3]) * 4.184
                        # Cv_corr might by NaN, this is a bug of gaussian
                        if math.isnan(Cv_corr):
                            Cv_corr = Cv
                        Cv_corr_list.append(Cv_corr)
                    else:
                        Cv_corr_list.append(Cv)
            return T_list, scale_list, Cv_list, Cv_corr_list

        import numpy as np

        if logs is None:
            logs = self.logs
        Cv_T = {}
        Cv_corr_T = {}
        for log in logs:
            results = process_log(log)
            if results is not None:
                T_list, scale_list, Cv_list, Cv_corr_list = results
                for i, T in enumerate(T_list):
                    if T not in Cv_T.keys():
                        Cv_T[T] = []
                    Cv_T[T].append(Cv_list[i])
                    if T not in Cv_corr_T.keys():
                        Cv_corr_T[T] = []
                    Cv_corr_T[T].append(Cv_corr_list[i])

        if len(Cv_T) == 0:
            raise Exception('All jobs are failed')

        Cv_ave = {}
        Cv_corr_ave = {}
        for T in Cv_T.keys():
            Cv_ave[T] = [np.mean(Cv_T[T]), np.std(Cv_T[T], ddof=1)]
            Cv_corr_ave[T] = [np.mean(Cv_corr_T[T]), np.std(Cv_corr_T[T], ddof=1)]

        return {
            'Cv'          : Cv_ave,
            'Cv-corrected': Cv_corr_ave,
        }

    def clean(self):
        pass
