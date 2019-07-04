import os, shutil
from .gmx import GmxSimulation
from ...wrapper.ppf import delta_ppf


class NvtExample(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'nvt-example'
        self.logs = ['nvt.log']
        self.n_atom_default = 1
        self.n_mol_default = 1

    def build(self, export=True, ppf=None):
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, 'init.pdb', length=28, silent=True)

        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd, mol_corr='init.pdb', length=30)
        if export:
            self.export(ppf=ppf)

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=298, jobname=None,
                dt=0.002, nst_run=int(1E6), nst_edr=10000, nst_trr=10000, nst_xtc=1000,
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

        # NVT production with Langevin thermostat
        self.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-nvt.mdp', T=T,
                                           dt=dt, nsteps=nst_run, nstenergy=nst_edr, nstxout=nst_trr, nstvout=nst_trr,
                                           nstxtcout=nst_xtc, restart=True)
        cmd = self.gmx.grompp(mdp='grompp-nvt.mdp', gro='conf.gro', top=top, tpr_out='nvt.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='nvt', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def analyze(self, dirs=None):
        pass
