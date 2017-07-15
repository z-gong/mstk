import os
import shutil

from ..simulation import Simulation
from ...errors import GmxError
from ...wrapper import GMX


class GmxSimulation(Simulation):
    def __init__(self, packmol_bin=None, dff_root=None, gmx_bin=None, jobmanager=None):
        super().__init__(packmol_bin=packmol_bin, dff_root=dff_root, jobmanager=jobmanager)
        self.gmx = GMX(gmx_bin=gmx_bin)
        self.logs = []  # used for checking whether the job is successfully finished

    def export(self, gro_out='conf.gro', top_out='topol.top', mdp_out='grompp.mdp',
               ppf=None, ff='TEAM_LS', minimize=False, vacuum=False):
        print('Generate GROMACS files ...')
        if ppf is not None:
            self.dff.typing([self.msd])  # in order to set the atom type
            self.dff.set_charge([self.msd], ppf)
            self.dff.export_gmx(self.msd, ppf, gro_out, top_out, mdp_out)
        else:
            self.dff.checkout([self.msd], table=ff)
            self.dff.export_gmx(self.msd, ff + '.ppf', gro_out, top_out, mdp_out)

        if minimize:
            print('Energy minimize ...')
            self.gmx.minimize(gro_out, top_out, name='em', silent=True, vacuum=vacuum, nprocs=self.jobmanager.nprocs)

            if os.path.exists('em.gro'):
                shutil.move('em.gro', gro_out)
            else:
                raise GmxError('Energy minimization failed')

    def check_finished(self, logs=None):
        if logs == None:
            logs = self.logs
        for log in logs:
            if not os.path.exists(log):
                return False
            with open(log) as f:
                lines = f.readlines()
            try:
                last_line = lines[-1]
            except:
                return False
            if not last_line.startswith('Finished mdrun'):
                return False
        return True
