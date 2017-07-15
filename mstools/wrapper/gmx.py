import os
import shutil
import subprocess
from subprocess import Popen, PIPE

from ..errors import GmxError


class GMX:
    TEMPLATE_DIR = os.path.abspath(os.path.dirname(__file__) + os.sep + '../template/gmx/')
    '''
    wrappers for GROMACS
    '''
    pass

    def __init__(self, gmx_bin):
        self.GMX_BIN = gmx_bin

    def grompp(self, mdp='grompp.mdp', gro='conf.gro', top='topol.top', tpr_out='md.tpr',
               cpt=None, maxwarn=3, silent=False, get_cmd=False):
        cmd = '%s -quiet -nobackup grompp -f %s -c %s -p %s -o %s -maxwarn %i' % (
            self.GMX_BIN, mdp, gro, top, tpr_out, maxwarn)
        if cpt != None:
            cmd = '%s -t %s' % (cmd, cpt)
        if get_cmd:
            return cmd
        else:
            (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
            sp = Popen(cmd.split(), stdout=stdout, stderr=stderr)
            sp.communicate()

    def mdrun(self, name='md', nprocs=1, rerun: str = None, extend=False, silent=False, get_cmd=False):
        cmd = '%s -quiet -nobackup mdrun -ntomp %i -deffnm %s' % (self.GMX_BIN, nprocs, name)

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

    def minimize(self, gro, top, nprocs=1, silent=False, name='em', vacuum=False):
        if not vacuum:
            self.prepare_mdp_from_template('t_em.mdp')
        else:
            self.prepare_mdp_from_template('t_em_vacuum.mdp')

        self.grompp(gro=gro, top=top, tpr_out=name + '.tpr', silent=silent)
        self.mdrun(name=name, nprocs=nprocs, silent=silent)

    def dos(self, trr, tpr, T, group='System', log_out='dos.log', get_cmd=False, silent=False):
        cmd = '%s -quiet -nobackup dos -f %s -s %s -T %f -g %s' % (self.GMX_BIN, trr, tpr, T, log_out)
        if get_cmd:
            cmd = 'echo "%s" | %s' % (group, cmd)
            return cmd
        else:
            (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
            sp = Popen(cmd.split(), stdin=PIPE, stdout=stdout, stderr=stderr)
            sp.communicate(input=group.encode())

    @staticmethod
    def prepare_mdp_from_template(template: str, mdp_out='grompp.mdp', T=298, P=1, nsteps=10000, dt=0.001,
                                  nstenergy=100, nstxout=0, nstvout=0, nstxtcout=10000, xtcgrps='System',
                                  restart=False, tcoupl='langevin', pcoupl='parrinello-rahman', constraints='h-bonds'):
        template = os.path.join(GMX.TEMPLATE_DIR, template)
        if not os.path.exists(template):
            raise GmxError('mdp template not found')

        if tcoupl.lower() == 'langevin':
            integrator = 'sd'
            tcoupl = 'no'
            tau_t = '2'
        elif tcoupl.lower() == 'nose-hoover':
            integrator = 'md-vv'
            tcoupl = 'nose-hoover'
            tau_t = '0.5'
        elif tcoupl.lower() == 'v-rescale':
            integrator = 'md'
            tcoupl = 'v-rescale'
            tau_t = '0.1'
        else:
            raise Exception('tcoupl not good, should be one of langvein, nose-hoover, v-rescale')

        if pcoupl.lower() == 'berendsen':
            tau_p = '1'
        elif pcoupl.lower() == 'parrinello-rahman':
            tau_p = '5'
        elif pcoupl.lower() == 'mttk':
            tau_p = '5'
            constraints = 'none'
        else:
            raise Exception('pcoupl not good, should be one of berendsen, parrinello-rahman, mttk')

        if restart:
            genvel = 'no'
            continuation = 'yes'
        else:
            genvel = 'yes'
            continuation = 'no'

        with open(template) as f_t:
            with open(mdp_out, 'w') as f_mdp:
                f_mdp.write(
                    f_t.read().replace('%T%', str(T)).replace('%P%', str(P)).replace('%nsteps%', str(int(nsteps))) \
                        .replace('%dt%', str(dt)).replace('%nstenergy%', str(nstenergy)) \
                        .replace('%nstxout%', str(nstxout)).replace('%nstvout%', str(nstvout)) \
                        .replace('%nstxtcout%', str(nstxtcout)).replace('%xtcgrps%', str(xtcgrps)) \
                        .replace('%genvel%', genvel).replace('%continuation%', continuation) \
                        .replace('%integrator%', integrator).replace('%tcoupl%', tcoupl).replace('%tau-t%', tau_t) \
                        .replace('%pcoupl%', pcoupl).replace('%tau-p%', tau_p) \
                        .replace('%constraints%', constraints))

    def energy(self, edr, properties: [str], begin=0, get_cmd=False):
        property_str = '\\n'.join(properties)
        cmd = '%s -quiet -nobackup energy -f %s -b %s' % (self.GMX_BIN, edr, str(begin))
        if get_cmd:
            cmd = 'echo "%s" | %s' % (property_str, cmd)
            return cmd
        else:
            sp = Popen(cmd.split(), stdout=PIPE, stdin=PIPE, stderr=PIPE)
            sp_out = sp.communicate(input=property_str.encode())[0]
            return sp_out

    def get_property(self, edr, property: str, begin=0) -> float:
        sp_out = self.energy(edr, properties=[property], begin=begin)

        for line in sp_out.decode().splitlines():
            if line.lower().startswith(property.lower()):
                return float(line.split()[1])
        raise GmxError('Invalid property')

    def get_property_and_stderr(self, edr, property: str, begin=0) -> [float, float]:
        sp_out = self.energy(edr, properties=[property], begin=begin)

        for line in sp_out.decode().splitlines():
            if line.lower().startswith(property.lower()):
                return (float(line.split()[1]), float(line.split()[2]))
        raise GmxError('Invalid property')


    def get_box(self, edr, begin=0) -> [float]:
        sp = subprocess.Popen([self.GMX_BIN, 'energy', '-f', edr, '-b', str(begin)], stdout=PIPE, stdin=PIPE,
                              stderr=PIPE)
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

    @staticmethod
    def scale_box(gro, gro_out, new_box: [float]):
        '''
        Scale gro box to desired size.
        The coordinate of all atoms are scaled.
        The velocities are not modified.
        Only support rectangular box
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
        shutil.copy(itp, 'topol.top')
        with open('topol.top', 'a') as f:
            f.write('\n[system]\n[molecules]\n')
            for i, molecule in enumerate(molecules):
                f.write('%s %i\n' % (molecule, numbers[i]))

    def pdb2gro(self, pdb, gro_out, box: [float]):
        if len(box) != 3:
            raise GmxError('Invalid box')
        sp = Popen([self.GMX_BIN, 'editconf', '-f', pdb, '-o', gro_out, '-box', str(box[0]), str(box[1]), str(box[2])])
        sp.communicate()

    def velacc(self, trr, tpr=None, group=None, begin=0, xvg_out='velacc', silent=False):
        if tpr == None:
            tpr = trr
        if group == None:
            raise GmxError('No group specifed')

        (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
        sp = Popen([self.GMX_BIN, 'velacc', '-f', trr, '-s', tpr, '-o', xvg_out, '-b', str(begin), '-mol',
                    '-nonormalize'],
                   stdin=PIPE, stdout=stdout, stderr=stderr)
        sp.communicate(input=str(group).encode())

    def msd_com(self, xtc, tpr, resname, beginfit=-1, endfit=-1, xvg_out=None, silent=False):
        ndx = 'com-' + resname + '.ndx'
        GMX.select_com(tpr, resname, ndx_out=ndx)
        if xvg_out == None:
            xvg_out = 'msd-com-%s.xvg' % resname

        (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
        sp = Popen([self.GMX_BIN, 'msd', '-f', xtc, '-s', tpr, '-n', ndx, '-o', xvg_out, '-nomw',
                    '-beginfit', str(beginfit), '-endfit', str(endfit)], stdout=stdout, stderr=stderr)
        sp_out = sp.communicate()[0]

        for line in sp_out.decode().splitlines():
            if line.startswith('D['):
                return line
        raise GmxError('Error running gmx msd')

    def traj_com(self, xtc, tpr, xtc_out=None, silent=False):
        ndx = 'com.ndx'
        self.select_com(tpr, 'all', ndx_out=ndx)
        if xtc_out == None:
            xtc_out = 'com.xtc'

        (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
        sp = Popen([self.GMX_BIN, 'traj', '-s', tpr, '-f', xtc, '-oxt', xtc_out, '-mol', '-n', ndx],
                   stdout=stdout, stderr=stderr)
        sp.communicate()
        return xtc_out, ndx

    def select_com(self, tpr, resname, ndx_out='selection.ndx'):
        sp = Popen([self.GMX_BIN, 'select', '-s', tpr, '-on', ndx_out], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        if resname == 'all':
            select_com_str = 'res_com of all'
        else:
            select_com_str = 'res_com of resname %s' % resname
        sp.communicate(input=select_com_str.encode())

    @staticmethod
    def generate_top_for_hvap(top, top_out):
        with open(top) as f:
            lines = f.read().splitlines()
        lines = [l for l in lines if not (l.startswith(';') or l == '')]

        f_out = open(top_out, 'w')

        line_number_molecule = []
        line_number_atom = []
        line_number_system = 0
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

        for n in range(line_number_system, len(lines)):
            f_out.write(lines[n] + '\n')

    def slice_gro_from_traj(self, trr, tpr, gro_out, begin, end, dt, silent=True):
        cmd = '%s -quiet -nobackup trjconv -f %s -s %s -o %s -b %s -e %s -dt %s -sep -pbc whole' % (
            self.GMX_BIN, trr, tpr, gro_out, str(begin), str(end), str(dt))
        (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
        sp = Popen(cmd.split(), stdin=PIPE, stdout=stdout, stderr=stderr)
        sp.communicate(input='System'.encode())

    def extend_tpr(self, tpr, extend, silent=True):
        cmd = '%s -quiet convert-tpr -s %s -o %s -extend %s' % (self.GMX_BIN, tpr, tpr, str(extend))
        (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
        sp = Popen(cmd.split(), stdout=stdout, stderr=stderr)
        sp.communicate()

    def get_length_of_traj(self, traj) -> float:
        cmd = '%s -quiet check -f %s' % (self.GMX_BIN, traj)
        sp = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
        out, err = sp.communicate()

        # The output of gmx check are in stderr
        for line in err.decode().splitlines():
            if line.startswith('Last frame'):
                length = float(line.split()[-1])
                return length

        raise Exception('Cannot open trajectory')

    @staticmethod
    def generate_gpu_multidir_cmds(dirs: [str], commands: [str], n_parallel=8, n_gpu=2) -> [[str]]:
        import math, re
        def replace_gpu_multidir_cmd(dirs: [str], cmd: str) -> str:
            if cmd.startswith('export'):
                pass
            elif cmd.find('gmx') != -1 and cmd.find(' mdrun ') != -1:
                cmd = re.sub('-ntomp\s+[0-9]+', '', cmd)  # remove -ntomp xx
                cmd = 'mpirun -np %i %s' % (len(dirs), cmd)  # add mpirun -np xx
                cmd += ' -multidir ' + ' '.join(dirs)  # add -multidir xx xx xx
                cmd += ' -gpu_id ' + ''.join(map(str, range(n_gpu))) * (len(dirs) // n_gpu) \
                       + ''.join(map(str, range(len(dirs) % n_gpu)))  # add -gpu_id 01230123012
            else:
                cmd = 'for i in %s; do cd $i; %s; done' % (' '.join(dirs), cmd)  # do it in every directory
            return cmd

        commands_list: [[str]] = []
        n_group: int = math.ceil(len(dirs) / n_parallel)
        for n in range(n_group):
            multidir = dirs[n * n_parallel:(n + 1) * n_parallel]
            commands_multidir: [str] = []
            for j, dirname in enumerate(multidir):
                commands_multidir.append('dir%i=%s' % (j, dirname))  # replace full dir path with $dir1, $dir2 ...
            for cmd in commands:
                commands_multidir.append(replace_gpu_multidir_cmd(['$dir%i' % j for j in range(len(multidir))], cmd))
            commands_list.append(commands_multidir)
        return commands_list
