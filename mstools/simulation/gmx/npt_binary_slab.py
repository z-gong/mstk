import os
import shutil

import numpy

from mstools.analyzer.series import is_converged
from .gmx import GmxSimulation
from ...unit import Unit


class NptBinarySlab(GmxSimulation):
    requirement = []


    def build(self):
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, 'init.pdb', length=self.length - 2, slab=True, silent=True)

        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd, mol_corr='init.pdb',
                                             length=self.length)

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=None, P=None, jobname=None, **kwargs):
        if os.path.abspath(model_dir) != os.getcwd():
            shutil.copy(os.path.join(model_dir, gro), gro)
            shutil.copy(os.path.join(model_dir, top), top)
            for f in os.listdir(model_dir):
                if f.endswith('.itp'):
                    shutil.copy(os.path.join(model_dir, f), '.')

        nprocs = self.jobmanager.nprocs
        commands = []
        self.gmx.prepare_mdp_from_template('t_npt.mdp', T=T, P=P / Unit.bar, nsteps=int(1E6))
        cmd = self.gmx.grompp(gro=gro, top=top, tpr_out=self.procedure, get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name=self.procedure, nprocs=nprocs, get_cmd=True)
        commands.append(cmd)
        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)

    def analyze(self):
        import panedr
        df = panedr.edr_to_df(self.procedure + '.edr')
        potential_series = df.Potential[500:]
        density_series = df.Density[500:]

        converged = is_converged(potential_series)
        if converged:
            converged = is_converged(density_series)

        return converged, {
            'potential': numpy.mean(potential_series),
            'density': numpy.mean(density_series)
        }
