import os
from mstk.utils import random_string


class Gauss:
    '''
    Wrapper for Gaussian for quantum calculations

    Parameters
    ----------
    gauss_bin : str
        Path of the Gaussian binary executable
    scrdir : str, optional
        Path of the scratch dir for Gaussian calculations.
        If not provided, the ramdisk `/dev/shm` will be used.
    '''

    _TEMPLATE_DIR = os.path.abspath(os.path.dirname(__file__) + os.sep + '../template/gauss/')

    def __init__(self, gauss_bin, scrdir=None):
        if scrdir is None:
            scrdir = '/dev/shm'
        self.GAUSS_BIN = gauss_bin
        self.GAUSS_ROOT = os.path.dirname(gauss_bin)
        self.SCRDIR = scrdir

    def generate_gjf_cv(self, name, py_mol, method='B3LYP', basis='6-31G*', scale=0.9613,
                        T_list=None):
        '''
        Genrate GJF input file for calculating heat capacity with Gaussian

        The molecule is provided through a pybel.Molecule object.
        Therefore, openbabel should be installed.
        The pybel.Molecule object should contain information about positions of all atoms.

        The hindered-rotor model is used to describing the low-frequency dihedral vibrations.

        A list of temperature can be provided.
        Then the thermodynamic analysis will be performed at all the provided temperatures,
        so that the result at intermediate temperature can be interpolated.

        Parameters
        ----------
        name : str
        py_mol : pybel.Molecule
        method : str
        basis : str
        scale : float
        T_list : list of float, optional
        '''
        if T_list is None:
            T_list = [298]
        gjf = name + '.gjf'
        with open(gjf, 'w') as f:
            f.write('%%chk=%(name)s.chk\n'
                    '# opt freq=hindrot %(method)s %(basis)s scale=%(scale).4f temperature=%(T).2f\n'
                    '\n'
                    'Title\n'
                    '\n'
                    '%(charge)i %(multiplicity)i\n'
                    % ({'name'        : name,
                        'method'      : method,
                        'basis'       : basis,
                        'scale'       : scale,
                        'T'           : T_list[0],
                        'charge'      : py_mol.charge,
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
                        % ({'name' : name,
                            'scale': scale,
                            'T'    : T
                            })
                        )

    def run_gjf(self, gjf, nprocs=1, memMB=1000, get_cmd=False, sh_out=None):
        '''
        Process the GJF file so that it is ready for multi-core Gaussian calculation, and then return the commands for calculation.

        The commands return can be feed into a :class:`~mstk.scheduler.Scheduler` for running on a HPC.

        Parameters
        ----------
        gjf : str
            The GJF file to be processed
        nprocs : str
            Number of cores to use
        memMB : int
            Maxium size of memory to use in unit of MB
        get_cmd : bool
            Whether or not return the commands for running Gaussian calculation
        sh_out : str, optional
            If provided, the commmands for running Gaussian calculation will be written to this shell script

        Returns
        -------
        commands : list of str, or None
            If `get_cmd` is True, the commands for running Gaussian will be returned
        '''
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
