import os
import shutil
import subprocess
from subprocess import Popen, PIPE

from ..errors import GaussError
from ..utils import create_mol_from_smiles, random_string


class Gauss:
    TEMPLATE_DIR = os.path.abspath(os.path.dirname(__file__) + os.sep + '../template/gauss/')
    '''
    wrappers for Gaussian
    '''
    pass

    def __init__(self, gauss_bin, scrdir=None):
        if scrdir is None:
            scrdir = '/dev/shm'
        self.GAUSS_BIN = gauss_bin
        self.GAUSS_ROOT = os.path.dirname(gauss_bin)
        self.SCRDIR = scrdir

    def generate_gjf_cv(self, name, py_mol, method='B3LYP', basis='6-31G*', scale=0.9613, T_list=[298]):
        gjf = name + '.gjf'
        with open(gjf, 'w') as f:
            f.write('%%chk=%(name)s.chk\n'
                    '# opt freq=hindrot %(method)s %(basis)s scale=%(scale).4f temperature=%(T).2f\n'
                    '\n'
                    'Title\n'
                    '\n'
                    '%(charge)i %(multiplicity)i\n'
                    % ({'name': name,
                        'method': method,
                        'basis': basis,
                        'scale': scale,
                        'T': T_list[0],
                        'charge': py_mol.charge,
                        'multiplicity': py_mol.spin
                        })
                    )
            for atom_line in py_mol.write('xyz').splitlines()[2:]:
                f.write(atom_line + '\n')
            f.write('\n')

            for T in T_list[1:]:
                f.write('--Link1--\n'
                        '%%chk=%(name)s.chk\n'
                        '# freq=(readfc,hindrot) geom=allcheck scale=%(scale).4f temperature=%(T).2f\n'
                        '\n'
                        % ({'name': name,
                            'scale': scale,
                            'T': T
                            })
                        )

    def run_gjf(self, gjf, nprocs=1, memMB: int = 1000, get_cmd=False, sh_out=None):
        with open(gjf) as f:
            content = f.read()
        with open(gjf, 'w') as f:
            f.write('%%nprocshared=%i\n' % nprocs)
            f.write('%%mem=%iMB\n' % memMB)
            f.write(content)

        JOB_DIR = os.path.join(self.SCRDIR, random_string(6))
        cmds = []
        cmds.append('JOB_DIR=' + JOB_DIR)
        cmds.append('mkdir -p ${JOB_DIR}')
        cmds.append('export GAUSS_EXEDIR=%s:%s/bsd' % (self.GAUSS_ROOT, self.GAUSS_ROOT))
        cmds.append('export GAUSS_SCRDIR=${JOB_DIR}')
        cmds.append('%s %s' % (self.GAUSS_BIN, gjf))
        cmds.append('rm -rf ${JOB_DIR}')

        if sh_out != None:
            with open(sh_out, 'w') as f:
                f.write('\n'.join(cmds))
        if get_cmd:
            return cmds
