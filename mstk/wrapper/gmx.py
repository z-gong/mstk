import os
import shutil
import subprocess
from subprocess import Popen, PIPE

from ..errors import GmxError
from ..utils import random_string


class GMX:
    '''
    Wrapper for GROAMCS for pre-processing, running and post-processing MD simulations.

    GROMACS version 2016 is recommended and has been tested thoroughly.
    Though this wrapper support both version 2016, 2018 and 2019.

    It's better to compile GROMACS in two binaries.
    One is the normal `gmx` binary without any acceleration.
    This is used for pre- and post-processing, which normally run on the front node of a HPC.
    The other one is the `mdrun` binary with full accelerations through MPI, OpenMP, GPU, etc...
    This will be submitted to :class:`~mstk.scheduler.Scheduler` to run on the compute node of a HPC.
    Refer to GROMACS documentation for the compiling options.

    If gmx_bin and gmx_mdrun are provided separately (recommended), you should make sure that they have the same version.

    Parameters
    ----------
    gmx_bin : str
        Path of the serial version gmx binary executable.
    gmx_mdrun : str, optional
        Path of the parallel/gpu version mdrun binary executable.
        If not provided, `gmx_bin` will be used also for mdrun.
    version : str, optional
        The version of GROMACS, e.g. 2016.6, 2019.6.
        If not specified, then gmx_bin will be called to check the version.
        In this case, make sure gmx_bin is available on this machine.

    Attributes
    ----------
    version : str
    majorversion : str
    '''

    _TEMPLATE_DIR = os.path.abspath(os.path.dirname(__file__) + os.sep + '../template/gmx/')

    def __init__(self, gmx_bin, gmx_mdrun=None, version=None):
        self.GMX_BIN = gmx_bin
        self.GMX_MDRUN = gmx_mdrun or gmx_bin + ' mdrun'
        if version is None:
            self._check_version()
        else:
            self.version = version
            self.majorversion = self.version.split('.')[0]

        # TODO temporary hack for dielectric constant in mdp
        self._DIELECTRIC = 1.0
        # TODO temporary hack for LJ96 function
        self._LJ96 = False

    def _check_version(self):
        cmd = '%s --version' % self.GMX_BIN
        try:
            sp = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
        except:
            raise GmxError('gmx not valid')
        stdout = sp.communicate()[0]
        for line in stdout.decode().splitlines():
            if line.startswith('GROMACS version'):
                break
        else:
            raise GmxError('gmx not valid')

        self.version = line.strip().split()[-1]
        if not (self.version.startswith('2016') or self.version.startswith(
                '2018') or self.version.startswith('2019')):
            raise GmxError('Supported GROMACS versions: 2016.x, 2018.x, 2019.x')

        self.majorversion = self.version.split('.')[0]

    def grompp(self, mdp='grompp.mdp', gro='conf.gro', top='topol.top', tpr_out='md.tpr',
               cpt=None, maxwarn=3, silent=False, get_cmd=False):
        '''
        Run `gmx grompp` directly or generate the command for `gmx grompp`.

        Refer to GROMACS documentation for details.

        Parameters
        ----------
        mdp : str
        gro : str
        top : str
        tpr_out : str
        cpt : str
        maxwarn : int
        silent : bool
            If set to True and `get_cmd` set to False, then `gmx grompp` will be executed silently, without output on screen.
        get_cmd : bool
            If set to True, the command for running `gmx grompp` will be returned, which can be feed into :class:`~mstk.scheduler.Scheduler`.
            If set to False, `gmx grompp` will be executed.

        Returns
        -------
        command : str or None
            If get_cmd set to True, will return the command for running `gmx grompp`.
            If get_cmd set to False, will return None.

        '''
        cmd = '%s -quiet -nobackup grompp -f %s -c %s -p %s -o %s -maxwarn %i' % (
            self.GMX_BIN, mdp, gro, top, tpr_out, maxwarn)
        if cpt is not None:
            cmd = '%s -t %s' % (cmd, cpt)
        if get_cmd:
            return cmd
        else:
            (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
            sp = Popen(cmd.split(), stdout=stdout, stderr=stderr)
            sp.communicate()

    def mdrun(self, name='md', nprocs=1, n_omp=None, rerun=None, extend=False, silent=False,
              get_cmd=False):
        '''
        Run `mdrun` directly or generate the command for `mdrun`.

        `mdrun` will use both MPI and OpenMP for hybrid parallelization.

        Refer to GROMACS documentation for details.

        Parameters
        ----------
        name : str
            Argument for `mdrun -deffnm`
        nprocs : int
            Number of CPU cores to use
        n_omp : int, optional
            Number of OpenMP threads to use.
            If not set, will be determined automatically from nprocs.
        rerun : str, optional
            If an XTC or TRR file privided, will perform `mdrun -rerun` on the trajectory.
        extend : bool
            If set to True, will perform `mdrun -cpi` on the checkpoint.
            A checkpoint named `name.cpt` should be avaiable in the current directory.
        silent : bool
            If set to True and `get_cmd` set to False, then `mdrun` will be executed silently, without output on screen.
        get_cmd : bool
            If set to True, the command for running `mdrun` will be returned, which can be feed into :class:`~mstk.scheduler.Scheduler`.
            If set to False, `mdrun` will be executed.

        Returns
        -------
        command : str or None
            If get_cmd set to True, will return the command for running `mdrun`.
            If get_cmd set to False, will return None.

        '''
        # TODO temporary hack for LJ96 function
        if self._LJ96:
            n_omp = 1
        ###

        if n_omp is None:  # n_omp is None means auto
            for i in [6, 4, 2]:  # available OpenMP threads: 6, 4, 2
                if nprocs % i == 0:
                    n_omp = i
                    break
            else:
                n_omp = nprocs
        n_mpi = nprocs // n_omp

        # TODO temporary hack for LJ96 function
        if self._LJ96:
            if name.find('hvap') != -1:
                n_mpi = 1
        ###

        cmd = '%s -quiet -nobackup -ntomp %i -deffnm %s' % (self.GMX_MDRUN, n_omp, name)
        # always use mpirun even if only one process
        cmd = 'mpirun -np %i ' % n_mpi + cmd

        if rerun is not None:
            cmd = '%s -rerun %s' % (cmd, rerun)

        if extend:
            cmd = '%s -cpi %s' % (cmd, name + '.cpt')

        if get_cmd:
            return cmd
        else:
            (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
            sp = Popen(cmd.split(), stdout=stdout, stderr=stderr)
            sp.communicate()

    def prepare_mdp_from_template(self, template, mdp_out='grompp.mdp', T=298, P=1, nsteps=10000,
                                  dt=0.001, TANNEAL=800,
                                  nstenergy=100, nstxout=0, nstvout=0, nstxtcout=10000,
                                  xtcgrps='System',
                                  restart=False, tcoupl='langevin', pcoupl='parrinello-rahman',
                                  constraints='h-bonds', ppm=0, dielectric=None):
        '''
        Generate MDP control script from template library for running GROMACS simulation.

        Refer to GROMACS documentation for the meaning of parameteres.

        The coupling time constant are determined automatically from `dt`, `tcoupl` and `pcoupl`.
        The neighbour list update frequency is determined automatically from `dt`.

        Parameters
        ----------
        template : str
            File name of the template to use.
            Take a look of the template directory to get available choices.
        mdp_out : str
        T : float
        P : float
        nsteps : int
        dt : float
        TANNEAL : float, optional
            The temperature for annealing.
            This argument should be used together with template `t_nvt_anneal.mdp`
            It is ignored for other templates.
        nstenergy : int
        nstxout : int
        nstvout : int
        nstxtcout : int
        xtcgrps : int
        restart : bool
            If set to True, grompp will not generated velocity.
        tcoupl : ['langevin', 'v-rescale', 'nose-hoover']
        pcoupl : ['parrinello-rahman', 'berendsen']
        constraints : ['h-bonds', 'none', 'all-bonds']
        ppm : float, optional
            The cosine acceleration strength for periodic perturbation simulation.
            This argument should be used together with template `t_npt_ppm.mdp`
            It is ignored for other templates.
        dielectric : float, optional

        '''
        template = os.path.join(GMX._TEMPLATE_DIR, template)
        if not os.path.exists(template):
            raise GmxError('mdp template not found')

        if tcoupl.lower() == 'langevin':
            integrator = 'sd'
            tcoupl = 'no'
            tau_t = str(0.001 / dt)  # inverse friction coefficient
        elif tcoupl.lower() == 'nose-hoover':
            integrator = 'md'
            tcoupl = 'nose-hoover'
            tau_t = '0.5'
        elif tcoupl.lower() == 'v-rescale':
            integrator = 'md'
            tcoupl = 'v-rescale'
            tau_t = '0.1'
        else:
            raise Exception('Invalid tcoupl, should be one of langvein, nose-hoover, v-rescale')

        if pcoupl.lower() == 'berendsen':
            tau_p = '1'
        elif pcoupl.lower() == 'parrinello-rahman':
            tau_p = '5'
        else:
            raise Exception('Invalid pcoupl, should be one of berendsen, parrinello-rahman')

        if restart:
            genvel = 'no'
            continuation = 'yes'
        else:
            genvel = 'yes'
            continuation = 'no'

        nstlist = max(1, int(0.01 / dt))

        if dielectric is None:
            dielectric = self._DIELECTRIC

        with open(template) as f_t:
            contents = f_t.read()
        contents = contents.replace('%T%', str(T)).replace('%P%', str(P)).replace('%nsteps%',
                                                                                  str(int(nsteps))) \
            .replace('%dt%', str(dt)).replace('%nstenergy%', str(nstenergy)) \
            .replace('%nstxout%', str(nstxout)).replace('%nstvout%', str(nstvout)) \
            .replace('%nstxtcout%', str(nstxtcout)).replace('%xtcgrps%', str(xtcgrps)) \
            .replace('%genvel%', genvel).replace('%continuation%', continuation) \
            .replace('%integrator%', integrator).replace('%tcoupl%', tcoupl).replace('%tau-t%',
                                                                                     tau_t) \
            .replace('%pcoupl%', pcoupl).replace('%tau-p%', tau_p) \
            .replace('%constraints%', constraints).replace('%TANNEAL%', str(TANNEAL)).replace(
            '%ppm%', str(ppm)) \
            .replace('%nstlist%', str(nstlist)).replace('%dielectric%', str(dielectric))

        with open(mdp_out, 'w') as f_mdp:
            f_mdp.write(contents)

    def energy(self, edr, properties, begin=0, end=None, fluct_props=False, get_cmd=False):
        '''
        Execute `gmx energy` program on a `edr` file.

        Refer to the GROMACS documentation for details.

        Parameters
        ----------
        edr : str
        properties : list of str
            The properties to be extracted.
            It corresponds to the input when calling gmx_energy manually.
            But instead of input one by one joined by Enter, the properties are provided in a list here.
        begin : float
        end : float, optional
        fluct_props : bool
        get_cmd : bool

        Returns
        -------
        out : str
            If set to True, the command for running `gmx energy` will be returned.
            If set to False, `gmx energy` will be executed silently and the output will be returned directly without cleanup.
        '''
        cmd = '%s -quiet -nobackup energy -f %s -b %s' % (self.GMX_BIN, edr, str(begin))
        if end is not None:
            cmd += ' -e %s' % (str(end))
        if fluct_props:
            cmd += ' -fluct_props'
        if get_cmd:
            property_str = '\\n'.join(properties)
            cmd = 'echo "%s" | %s' % (property_str, cmd)
            return cmd
        else:
            sp = Popen(cmd.split(), stdout=PIPE, stdin=PIPE, stderr=PIPE)
            property_str = '\n'.join(properties)
            out, err = sp.communicate(input=property_str.encode())
            return out.decode()

    def get_fluct_props(self, edr, begin=0, end=None):
        '''
        Get thermal expansion and compressibility from a `edr` file using fluctuation of enthalpy, volume.

        Parameters
        ----------
        edr : str
        begin : float
        end : float, optional

        Returns
        -------
        expansion : float
        compressibility : float
        '''
        sp_out = self.energy(edr, properties=['temp', 'vol', 'enthalpy'], begin=begin, end=end,
                             fluct_props=True)

        expansion = None
        compressibility = None
        for line in sp_out.splitlines():
            if line.startswith('Coefficient of Thermal Expansion Alpha_P'):
                expansion = float(line.split()[-2])
            elif line.startswith('Isothermal Compressibility Kappa'):
                compressibility = float(line.split()[-2]) * 1e5
        return expansion, compressibility

    def get_properties_stderr(self, edr, properties, begin=0, end=None):
        '''
        Extract properties from a `edr` file by calling `gmx energy` program.
        The difference between this method and :func:`energy` is that,
        :func:`energy` will return the raw output string without any cleanup,
        while this method will extract the values and standard errors from the output and return only the results.

        Parameters
        ----------
        edr : str
        properties : list of str
        begin : float
        end : float, optional

        Returns
        -------
        results : list of tuple of float
            Each element of this list is a tuple of two floats, which are average value and standard error of each property.
        '''
        sp_out = self.energy(edr, properties=properties, begin=begin, end=end)

        lines = sp_out.splitlines()
        results = []
        for prop in properties:
            for line in lines:
                if line.lower().startswith(prop.lower()):
                    results.append((float(line.split()[1]), float(line.split()[2])))
                    break
            else:
                raise GmxError('Invalid property')
        return results

    def get_property_stderr(self, edr, prop, begin=0, end=None):
        '''
        Extract the average value and standard error of a property from a `edr` file
        by calling `gmx energy` program.

        Parameters
        ----------
        edr : str
        prop: str
        begin : float
        end : float, optional

        Returns
        -------
        average : float
        stderr : float
        '''
        return self.get_properties_stderr(edr, [prop], begin, end)[0]

    def get_box(self, edr, begin=0):
        '''
        Extract the average box size from a `edr` file by calling `gmx energy` program.

        Only rectangular box is supported.

        Parameters
        ----------
        edr : str
        begin : float

        Returns
        -------
        lx : float
        ly : float
        lz : float
        '''
        sp = subprocess.Popen([self.GMX_BIN, 'energy', '-f', edr, '-b', str(begin)],
                              stdout=PIPE, stdin=PIPE, stderr=PIPE)
        sp_out = sp.communicate(input='Box-'.encode())[0]

        box = [0, 0, 0]
        for line in sp_out.decode().splitlines():
            if line.startswith('Box-X'):
                box[0] = float(line.split()[1])
            if line.startswith('Box-Y'):
                box[1] = float(line.split()[1])
            if line.startswith('Box-Z'):
                box[2] = float(line.split()[1])
        return box

    def density(self, trj, tpr, xvg='density.xvg', group='System', begin=0, end=None,
                center=False, silent=False, get_cmd=False):
        '''
        Execute `gmx density` to calculate the density profile in `z` direction on a trajectory.

        Refer to GROMACS documentation for details.

        Parameters
        ----------
        trj : str
        tpr : str
        xvg : str
        group : str
        begin : float
        end : float, optional
        center : bool
        silent : bool
            If set to True and `get_cmd` set to False, then `gmx density` will be executed silently, without output on screen.
        get_cmd : bool
            If set to True, the command for running `gmx density` will be returned, which can be feed into :class:`~mstk.scheduler.Scheduler`.
            If set to False, `gmx density` will be executed.

        Returns
        -------
        command : str or None
            If get_cmd set to True, will return the command for running `gmx_density`.
            If get_cmd set to False, will return None.
        '''
        cmd = '%s -quiet -nobackup density -f %s -s %s -b %f -o %s' % (
            self.GMX_BIN, trj, tpr, begin, xvg)
        if end is not None:
            cmd += ' -e %f' % end
        inp = group
        if center:
            cmd += ' -center'
            inp = '%s\n%s' % (group, group)

        if get_cmd:
            cmd = 'echo "%s" | %s' % (inp, cmd)
            return cmd
        else:
            (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
            sp = Popen(cmd.split(), stdin=PIPE, stdout=stdout, stderr=stderr)
            sp.communicate(input=inp.encode())

    @staticmethod
    def scale_box(gro, gro_out, new_box):
        '''
        Scale GRO box to desired size.

        The coordinate of all atoms are scaled, while the velocities are not modified.
        Only rectangular box is supported.

        Parameters
        ----------
        gro : str
        gro_out : str
        new_box : list of float
            The lengths of the desired rectangular box.
        '''
        NAtoms = 0
        nresidue = []
        residue = []
        element = []
        natom = []
        xyz = []
        vxyz = []
        box = []
        if not os.path.exists(gro):
            raise GmxError('gro not found')
        with open(gro) as f_gro:
            lines = f_gro.read().splitlines()

        for i, line in enumerate(lines):
            n = i + 1
            if n == 1:
                title = line
            if n == 2:
                NAtoms = int(line.strip())
                continue
            if n > 2 and n <= 2 + NAtoms:
                nresidue.append(int(line[:5]))
                residue.append(line[5:10])
                element.append(line[10:15])
                natom.append(int(line[15:20]))
                x = float(line[20:28])
                y = float(line[28:36])
                z = float(line[36:44])
                vx = float(line[44:52])
                vy = float(line[52:60])
                vz = float(line[60:68])
                xyz.append([x, y, z])
                vxyz.append([vx, vy, vz])
            if n == 3 + NAtoms:
                box = [float(word) for word in line.strip().split()[:3]]
                break

        scale = [new_box[i] / box[i] for i in range(3)]
        xyz = [[i[0] * scale[0], i[1] * scale[1], i[2] * scale[2]] for i in xyz]

        with open(gro_out, 'w') as f_out:
            f_out.write('Scaled : %s\n%i\n' % (title, NAtoms))
            for i in range(NAtoms):
                f_out.write('%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n'
                            % (nresidue[i], residue[i], element[i], natom[i],
                               xyz[i][0], xyz[i][1], xyz[i][2], vxyz[i][0], vxyz[i][1], vxyz[i][2]))
            f_out.write('%f %f %f\n' % (new_box[0], new_box[1], new_box[2]))

    @staticmethod
    def generate_top(itp, molecules, numbers):
        '''
        Generate a `topol` file that includes a existent `itp` file.

        Parameters
        ----------
        itp : str
        molecules : list of str
            Name of molecules in the system
        numbers : list of int
            Number of these molecules in the system
        '''
        shutil.copy(itp, 'topol.top')
        with open('topol.top', 'a') as f:
            f.write('\n[system]\n[molecules]\n')
            for i, molecule in enumerate(molecules):
                f.write('%s %i\n' % (molecule, numbers[i]))

    @staticmethod
    def get_top_mol_numbers(top):
        '''
        Get the name and number of molecules in a topol file from the `[ molecules ]` section.

        Parameters
        ----------
        top : str

        Returns
        -------
        mol_numbers : list of tuple
            Each element is a tuple containing a str and a int, which are the name and number of each molecule specie.

        '''
        with open(top) as f:
            lines = f.read().splitlines()

        newlines = []
        mols = []
        START = False
        for line in lines:
            if line.find('[') != -1 and line.find('molecules') != -1:  # [ molecules ]
                START = True
                newlines.append(line)
                continue
            if not START:
                newlines.append(line)
            if START and not line.strip() == '' and not line.startswith(';'):
                words = line.strip().split()
                mols.append((words[0], int(words[1])))
        return mols

    @staticmethod
    def modify_top_mol_numbers(top, numbers):
        '''
        Modify the number of molecules in a `topol` file at the `[ molecules ]` section.

        The length of numbers must be equal to the original length of numbers in the `topol` file.
        Otherwise, an Exception will be raised.

        Parameters
        ----------
        top : str
        numbers : list of int
        '''
        with open(top) as f:
            lines = f.read().splitlines()

        newlines = []
        mols = []
        START = False
        for line in lines:
            if line.find('[') != -1 and line.find('molecules') != -1:  # [ molecules ]
                START = True
                newlines.append(line)
                continue
            if not START:
                newlines.append(line)
            if START and not line.strip() == '' and not line.startswith(';'):
                words = line.strip().split()
                mols.append((words[0], int(words[1])))
        if len(mols) != len(numbers):
            raise GmxError('Type of molecules in top not consistent')

        for i, (mol_name, _) in enumerate(mols):
            newlines.append('%s\t%i' % (mol_name, numbers[i]))

        with open(top, 'w')  as f:
            f.write('\n'.join(newlines))

    def replicate_gro(self, gro, top, nbox, silent=True):
        '''
        Enlarge a GRO box by replicate it with `gmx genconf`.

        The numbers of molecules in the `topol` file will be updated accordingly.

        Parameters
        ----------
        gro : str
        top : str
        nbox : list of int
            Three element denotes the number of replicates in x, y  and z directions.
        silent : bool
        '''
        from functools import reduce
        import operator

        gro_tmp = random_string(8) + '.gro'
        cmd = '%s -quiet genconf -f %s -o %s -nbox %i %i %i' % (
            self.GMX_BIN, gro, gro_tmp, nbox[0], nbox[1], nbox[2])
        (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
        sp = Popen(cmd.split(), stdout=stdout, stderr=stderr)
        sp.communicate()
        shutil.move(gro_tmp, gro)

        with open(top) as f:
            content = f.read()

        LASTLINE = '\n'
        LAST = False
        for line in content.splitlines():
            if line.find('molecules') > -1 and line.find('[') > -1:
                LAST = True
                continue
            if LAST and not line.strip() == '' and not line.startswith(';'):
                LASTLINE += line + '\n'
        LASTLINE *= (reduce(operator.mul, nbox, 1) - 1)

        with open(top, 'a') as f:
            f.write(LASTLINE)

    def pdb2gro(self, pdb, gro_out, box, silent=False):
        '''
        Convert a PDB file to GRO file by calling `gmx editconf`.

        Parameters
        ----------
        pdb : str
        gro_out : str
        box : list of float
            The size of rectangular box.
        silent : bool
        '''
        if len(box) != 3:
            raise GmxError('Invalid box')

        (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
        sp = Popen(
            [self.GMX_BIN, 'editconf', '-f', pdb, '-o', gro_out, '-box', str(box[0]), str(box[1]),
             str(box[2])],
            stdin=PIPE, stdout=stdout, stderr=stderr)
        sp.communicate()

    def velacc(self, trr, tpr=None, group=None, begin=0, xvg_out='velacc', silent=False):
        '''
        Calculte the velocity auto correlation of atoms in a trajectory by calling `gmx velacc`.

        Refer to GROMACS documentation for details.

        Parameters
        ----------
        trr : str
        tpr : str
        group : str
        begin : float
        xvg_out : str
        silent : bool
        '''
        if tpr is None:
            tpr = trr
        if group is None:
            raise GmxError('No group specifed')

        (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
        sp = Popen(
            [self.GMX_BIN, 'velacc', '-f', trr, '-s', tpr, '-o', xvg_out, '-b', str(begin), '-mol',
             '-nonormalize'],
            stdin=PIPE, stdout=stdout, stderr=stderr)
        sp.communicate(input=str(group).encode())

    def diffusion(self, xtc, tpr, group='System', mol=False, begin=0, end=None, xvg_out='msd.xvg',
                  beginfit=-1, endfit=-1):
        '''
        Calculate the diffusion coefficient from a trajectory file by calling `gmx msd`.

        Refer to GROMACS documentation for details.

        Parameters
        ----------
        xtc : str
        tpr : str
        group : str
        mol : bool
        begin : float
        end : float, optional
        xvg_out : str
        beginfit : float, optional
        endfit : str, optional

        Returns
        -------
        diffusion : float
            Diffusion coefficient in unit of cm^2/s.
        stderr : float
            Standard error.
        '''
        cmd = '%s -quiet -nobackup msd -f %s -s %s -o %s -b %s -beginfit %s -endfit %s' % (
            self.GMX_BIN, xtc, tpr, xvg_out, str(begin), str(beginfit), str(endfit))
        if end is not None:
            cmd += ' -e %s' % str(end)
        if mol:
            # calculate the MSD of COM of molecules
            cmd += ' -mol'

        sp = Popen(cmd.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out, err = sp.communicate(input=str(group).encode())

        for line in out.decode().splitlines():
            if line.startswith('D['):
                words = line.strip().split(']')[-1].strip().split()
                diffusion = float(words[0])
                stderr = float(words[2][:-1])
                unit = float(words[3])
                diffusion *= unit  # cm^2/s
                stderr *= unit  # cm^2/s
                return diffusion, stderr

        raise GmxError('Error running gmx msd')

    def traj_com(self, xtc, tpr, trj_out='com.xtc', begin=0, end=0, silent=False):
        '''
        Generate a trajectory of the centers of mass of molecules in a trajectory.

        It calls :func:`select_com` to create a index file of the centers of mass of all molecules.
        Then `gmx traj` is called to calculate the positions of each center of mass and output the trajectory.

        Refer to the GROMACS documentation for details.

        Parameters
        ----------
        xtc : str
        tpr : str
        trj_out : str
        begin : float
        end : float
        silent : bool
        '''
        ndx = 'com.ndx'
        self.select_com(tpr, 'all', ndx_out=ndx)

        (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
        cmd = '%s -quiet -nobackup traj -s %s -f %s -oxt %s -com -mol -n %s -b %f -e %f' % (
            self.GMX_BIN, tpr, xtc, trj_out, ndx, begin, end)
        sp = Popen(cmd.split(), stdout=stdout, stderr=stderr)
        sp.communicate()
        return trj_out, ndx

    def select_com(self, tpr, resname='all', ndx_out='com.ndx'):
        '''
        Create an index file of the centers of mass of each molecule by calling `gmx select`.

        Refer to GROAMCS documentation for details.

        Parameters
        ----------
        tpr : str
        resname : str
        ndx_out : str
        '''
        cmd = '%s -quiet -nobackup select -s %s -on %s' % (self.GMX_BIN, tpr, ndx_out)
        sp = Popen(cmd.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
        if resname == 'all':
            select_com_str = 'res_com of all'
        else:
            select_com_str = 'res_com of resname %s' % resname
        sp.communicate(input=select_com_str.encode())

    @staticmethod
    def generate_top_for_hvap(top, top_out):
        '''
        Generate a new `topol` file for calculating the heat of vaporization from existent `topol` file for bulk simulation.

        Essentially, this will add a `[ exclusion ]` section in the `topol` file for each molecule.
        All intra-molecular interactions will be excluded.
        With this modified `topol` file, `gmx mdrun -rerun` on the bulk trajectory will calculate
        only the inter-molecular energies, which is approximately the cohesive energy.

        Parameters
        ----------
        top : str
        top_out : str
        '''

        with open(top) as f:
            lines = f.read().splitlines()
        lines = [l for l in lines if not (l.startswith(';') or l == '')]

        f_out = open(top_out, 'w')

        line_number_molecule = []
        line_number_atom = []
        line_number_system = None
        for n, line in enumerate(lines):
            if line.find('[') != -1 and line.find('moleculetype') != -1:
                line_number_molecule.append(n)
            if line.find('[') != -1 and line.find('atoms') != -1:
                line_number_atom.append(n)
            if line.find('[') != -1 and line.find('system') != -1:
                line_number_system = n

        n_molecules = len(line_number_molecule)

        for n in range(line_number_molecule[0]):
            f_out.write(lines[n] + '\n')

        for i in range(n_molecules):
            for n in range(line_number_molecule[i], line_number_atom[i]):
                f_out.write(lines[n] + '\n')
            line_number_next_section = line_number_molecule[i + 1] if i < n_molecules - 1 else line_number_system
            if line_number_next_section is None:
                line_number_next_section = len(lines)
            n_atoms = 0
            f_out.write('[ atoms ]\n')
            for n in range(line_number_atom[i] + 1, line_number_next_section):
                line = lines[n]
                if line.find('[') != -1 or line.startswith('#'):
                    f_out.write('[ bonds ]\n[ exclusions ]\n')
                    for i in range(1, n_atoms):
                        exclusions = range(i, n_atoms + 1)
                        f_out.write(' '.join(list(map(str, exclusions))) + '\n')
                    break
                f_out.write(line + '\n')
                n_atoms += 1

        if line_number_system is not None:
            for n in range(line_number_system, len(lines)):
                f_out.write(lines[n] + '\n')

    def extend_tpr(self, tpr, extend, silent=False):
        '''
        Modify the `nsteps` record in a `tpr` file by calling `gmx convert-tpr` for extending simulation.

        Refer to GROMACS documentation for details.

        Parameters
        ----------
        tpr : str
        extend : float
        silent : bool
        '''
        cmd = '%s -quiet convert-tpr -s %s -o %s -extend %s' % (self.GMX_BIN, tpr, tpr, str(extend))
        (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
        sp = Popen(cmd.split(), stdout=stdout, stderr=stderr)
        sp.communicate()

    def generate_gpu_multidir_cmds(self, dirs, commands, n_parallel, n_gpu=0, n_omp=None,
                                   n_procs=None):
        '''
        Generate commands for performing `multidir` simulation with GROMACS.

        multidir simulation enables several similar simulations being run at the same time in different directories.
        This is a powerful feature allowing for much better utilization of GPU resources.

        Parameters
        ----------
        dirs : list of str
        commands : list of str
        n_parallel : int
        n_gpu : int, optional
        n_omp : int, optional
            Set n_omp in most case.
            When n_procs is set, n_omp will have no effect.
        n_procs : int, optional

        Returns
        -------
        commands_list : list of list of str

        '''
        import math, re
        def replace_gpu_multidir_cmd(dirs: [str], cmd: str) -> str:
            n_multi = len(dirs)
            if cmd.startswith('export '):
                pass

            elif cmd.find('gmx') != -1 and cmd.find(' mdrun ') != -1:
                n_mpi = n_multi
                n_thread = n_omp

                if n_procs is not None:
                    if cmd.find('mpirun ') != -1 and cmd.find('mpirun -np 1 ') == -1:
                        # optimize n_mpi and n_thread
                        # GPU tasks prefer more n_mpi
                        # ensure n_mpi equals 2*n
                        # do not optimize for mpirun -np 1. This happens for rerun hvap
                        optimal_omp = [6, 4, 2] if n_gpu > 0 else [8, 6, 4, 2]
                        for i in optimal_omp:
                            if n_procs % (n_multi * i * 2) == 0:
                                n_thread = i
                                n_mpi = n_procs // n_thread
                                break

                    else:
                        # set n_mpi equal to n_multi
                        n_thread = n_procs // n_mpi

                cmd = re.sub(r'mpirun\s+-np\s+[0-9]+', '', cmd)  # remove mpirun -np xx
                cmd = re.sub(r'-ntomp\s+[0-9]+', '', cmd)  # remove -ntomp xx

                cmd = 'mpirun -np %i %s' % (n_mpi, cmd)  # add mpirun -np xx
                cmd += ' -multidir ' + ' '.join(dirs)  # add -multidir xx xx xx

                ### meaning of -gpu_id is changed in GROMACS 2018. Disable gpu assignment
                if n_gpu > 0 and self.majorversion == '2016':
                    cmd += ' -gpu_id ' + ''.join(map(str, range(n_gpu))) * (n_mpi // n_gpu) \
                           + ''.join(map(str, range(n_mpi % n_gpu)))  # add -gpu_id 01230123012
                if n_thread is not None:
                    cmd += ' -ntomp %i' % n_thread

            else:
                cmd = 'for i in %s; do cd $i; %s; done' % (
                    ' '.join(dirs), cmd)  # do it in every directory
            return cmd

        commands_list: [[str]] = []
        n_group: int = math.ceil(len(dirs) / n_parallel)
        for n in range(n_group):
            multidir = dirs[n * n_parallel:(n + 1) * n_parallel]
            commands_multidir: [str] = []
            for j, dirname in enumerate(multidir):
                commands_multidir.append(
                    'dir%i=%s' % (j, dirname))  # replace full dir path with $dir1, $dir2 ...
            for cmd in commands:
                commands_multidir.append(
                    replace_gpu_multidir_cmd(['$dir%i' % j for j in range(len(multidir))], cmd))
            commands_list.append(commands_multidir)
        return commands_list

    @staticmethod
    def modify_lj96(itp_files, top_files, mdp_files, xvg_files):
        '''
        Convert the current simulation files for LJ-12-6 vdW form into LJ-9-6 form.

        This is a trick to simplify the exporting of LJ-9-6 vdW form in GROMACS.

        Parameters
        ----------
        itp_files : None or list of str
            The `itp` files to be modified.
            If not set, all the `itp` files in the current directory will be modified.
        top_files : None or list of str
            The `top` files to be modified.
            If not set, all the `top` files in the current directory will be modified.
        mdp_files : None or list of str
            The `mdp` files to be modified.
            If not set, all the `mdp` files in the current directory will be modified.
        xvg_files : list of str
            The `xvg` files to be generated. They contain the tabulated LJ-9-6 potential and force.
        '''
        if itp_files is None:
            itp_files = list(filter(lambda x: x.endswith('.itp'), os.listdir('.')))
        for itp in itp_files:
            GMX._modify_lj96_itp(itp)

        if top_files is None:
            top_files = list(filter(lambda x: x.endswith('.top'), os.listdir('.')))
        for top in top_files:
            GMX._modify_lj96_top(top)

        if mdp_files is None:
            mdp_files = list(filter(lambda x: x.endswith('.mdp'), os.listdir('.')))
        for mdp in mdp_files:
            GMX._modify_lj96_mdp(mdp)

        for xvg in xvg_files:
            shutil.copy(os.path.join(GMX._TEMPLATE_DIR, 'table6-9.xvg'), xvg)

    @staticmethod
    def _modify_lj96_itp(itp):
        with open(itp) as f:
            lines = f.read().splitlines()
        new_lines = []
        for line in lines:
            if line.strip().split()[:3] == ['1', '2', 'yes']:
                new_lines.append(line + ' 9')
            else:
                new_lines.append(line)
        with open(itp, 'w') as f:
            f.write('\n'.join(new_lines))

    @staticmethod
    def _modify_lj96_top(top):
        with open(top) as f:
            lines = f.read().splitlines()

        new_lines = []
        ATOMLINE = False
        for line in lines:
            if line.startswith(';') or line.strip() == '':
                new_lines.append(line)
                continue

            if line.find('[') != -1 and line.find('atoms') != -1:
                ATOMLINE = True
                new_lines.append(line)
                continue

            if ATOMLINE and line.find('[') != -1:
                ATOMLINE = False

            if ATOMLINE:
                words = line.strip().split()
                words[5] = words[0]
                tmp = ''
                for word in words:
                    tmp += ' %10s' % word
                new_lines.append(tmp)
            else:
                new_lines.append(line)
                continue

        with open(top, 'w') as f:
            f.write('\n'.join(new_lines))

    @staticmethod
    def _modify_lj96_mdp(mdp):
        with open(mdp) as f:
            lines = f.read().splitlines()
        new_lines = []
        for line in lines:
            if line.strip().startswith('cutoff-scheme') or \
                    line.strip().startswith('cutoff_scheme') or \
                    line.strip().startswith('cutoffscheme') or \
                    line.strip().startswith('rlist'):
                pass
            else:
                new_lines.append(line)

            if line.strip().startswith('rvdw'):
                rvdw = line.split('=')[-1].strip()

        new_lines.append('\n; vdw table')
        new_lines.append('cutoff-scheme = group')
        new_lines.append('rlist = ' + rvdw)
        new_lines.append('vdwtype = user')
        with open(mdp, 'w') as f:
            f.write('\n'.join(new_lines))
