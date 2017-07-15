import os
import shutil

from .gmx import GmxSimulation


class NvtSlab(GmxSimulation):
    requirement = []

    def build(self, ppf=None, minimize=False):
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, 'init.pdb', length=self.length - 2, silent=True)

        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd, mol_corr='init.pdb',
                                         size=[self.length, self.length, self.length * 5])
        self.export(ppf=ppf, minimize=minimize)

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=None, jobname=None,
                dt=0.001, nst_eq=int(4E5), nst_run=int(1E6), nst_edr=100, nst_trr=int(1E4), nst_xtc=int(1E3), **kwargs):
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

        # NVT equilibrium
        self.gmx.prepare_mdp_from_template('t_nvt_slab.mdp', mdp_out='grompp-eq.mdp', T=T,
                                           dt=0.001, nsteps=nst_eq, nstxtcout=0)
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

    def analyze(self):
        import panedr
        df = panedr.edr_to_df(self.procedure + '.edr')
        surface_tension = df['#Surf*SurfTen']
        return {
            'surface_tension': surface_tension
        }
