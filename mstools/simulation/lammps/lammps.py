import math
import os
import shutil

from ..simulation import Simulation
from ...errors import LammpsError
from ...jobmanager import Local
from ...utils import create_mol_from_smiles
from ...wrapper import Packmol, DFF, Lammps
from ...unit import Unit


class LammpsSimulation(Simulation):
    def __init__(self, packmol_bin=None, dff_root=None, lmp_bin=None, jobmanager=Local()):
        super().__init__(packmol_bin=packmol_bin, dff_root=dff_root, jobmanager=jobmanager)
        self.LMP_BIN = lmp_bin

    def build(self, minimize=False):
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, 'init.pdb', length=self.length - 2, silent=True)

        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd, mol_corr='init.pdb',
                                         length=self.length)
        self.export(minimize=minimize)

    def export(self, ff='TEAM_LS', data_out='data', in_lmp='em.lmp', minimize=False):
        print('Export lammps files...')
        self.dff.checkout(['init.msd'], table=ff)
        self.dff.export_lammps('init.msd', ff + '.ppf', data_out, in_lmp)

        if minimize:
            lammps = Lammps(self.LMP_BIN)
            print('Energy minimizing...')
            lammps.run(in_lmp, silent=True)

            if os.path.exists('em.data'):
                shutil.move('em.data', data_out)
            else:
                raise LammpsError('Energy minimization failed')

    def prepare(self, model_dir='.', T=None, P=None):
        shutil.copy('build/em.data', 'em.data')
        special, bond, angle, dihedral, improper = Lammps.get_intra_style_from_lmp('build/em.lmp')
        Lammps.prepare_lmp_from_template('t_npt.lmp', 'in.lmp', 'em.data', T, P / Unit.bar, int(1E3),
                                         self.n_mol_list[0], special, bond, angle, dihedral, improper)
