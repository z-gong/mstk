import os

from ..simulation import Simulation
from ...wrapper import GMX
from ...utils import get_last_line


class GmxSimulation(Simulation):
    def __init__(self, packmol_bin=None, dff_root=None, gmx_bin=None, gmx_mdrun=None, jobmanager=None, **kwargs):
        super().__init__(packmol_bin=packmol_bin, dff_root=dff_root, jobmanager=jobmanager, **kwargs)
        self.gmx = GMX(gmx_bin=gmx_bin, gmx_mdrun=gmx_mdrun)
        self.logs = []  # used for checking whether the job is successfully finished

    def export(self, gro_out='conf.gro', top_out='topol.top', mdp_out='grompp.mdp', ppf=None, ff=None):
        print('Generate GROMACS files ...')
        msd = self.msd
        self.dff.set_formal_charge([msd])
        if ppf is not None:
            self.dff.typing([msd])  # in order to set the atom type
            self.dff.set_charge([msd], ppf)
            self.dff.export_gmx(msd, ppf, gro_out, top_out, mdp_out)
        else:
            ppf_out = 'ff.ppf'
            self.dff.checkout([msd], table=ff, ppf_out=ppf_out)
            self.dff.export_gmx(msd, ppf_out, gro_out, top_out, mdp_out)

    def fast_export_single(self, gro_out='conf.gro', top_out='topol.top', mdp_out='grompp.mdp', ppf=None, ff=None):
        print('Generate GROMACS files ...')
        msd = self._single_msd
        self.dff.set_formal_charge([msd])
        if ppf is not None:
            self.dff.typing([msd])  # in order to set the atom type
            self.dff.set_charge([msd], ppf)
            self.dff.export_gmx(msd, ppf, gro_out, top_out, mdp_out)
        else:
            ppf_out = 'ff.ppf'
            self.dff.checkout([msd], table=ff, ppf_out=ppf_out)
            self.dff.export_gmx(msd, ppf_out, gro_out, top_out, mdp_out)

    def check_finished(self, logs=None):
        if logs is None:
            logs = self.logs
        for log in logs:
            if not os.path.exists(log):
                return False
            try:
                last_line = get_last_line(log)
            except:
                return False
            if not last_line.startswith('Finished mdrun'):
                return False
        return True
