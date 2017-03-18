import os
import subprocess
from subprocess import PIPE, STDOUT
import shutil

class GMX:
    TEMPLATE_DIR = os.path.abspath(os.path.dirname(__file__) + os.sep + '../template/gmx/')
    '''
    wrappers for GROMACS
    '''
    pass

    @staticmethod
    def grompp(mdp='grompp.mdp', gro='conf.gro', top='topol.top', out='md.tpr', maxwarn=3):
        sp = subprocess.run(['gmx', 'grompp', '-f', mdp, '-c', gro, '-p', top, '-o', out, '-maxwarn', str(maxwarn)])

    @staticmethod
    def mdrun(tpr, nprocs=1):
        if nprocs == 1:
            sp = subprocess.run(['gmx', 'mdrun', '-v', '-deffnm', tpr])
        else:
            sp = subprocess.run(['mpirun', '-np', str(nprocs), 'gmx', 'mdrun', '-v', '-deffnm', tpr])


    @staticmethod
    def prepare_mdp_from_template(template: str, mdp='grompp.mdp', T=None, P=None, nsteps=None,
            nstenergy=1000, nstvout = 0,
            nstxtcout=10000, xtcgrps = 'System',
            genvel=True):

        if T == None or P == None or nsteps == None:
            raise Exception('T, P, nsteps are required')

        genvel = 'yes' if genvel else 'no'

        template = os.path.join(GMX.TEMPLATE_DIR, template)
        if not os.path.exists(template):
            raise Exception('mdp template not found')

        with open(template) as f_t:
            with open(mdp, 'w') as f_mdp:
                f_mdp.write(f_t.read().replace('%T%', str(T)).replace('%P%', str(P)).replace('%nsteps%', str(nsteps))\
                        .replace('%nstenergy%', str(nstenergy)).replace('%nstvout%', str(nstvout))\
                        .replace('%nstxtcout%', str(nstxtcout)).replace('%xtcgrps%', str(xtcgrps))\
                        .replace('%genvel%', str(genvel)))

    @staticmethod
    def get_property(edr, property_str, begin=0) -> float:
        sp = subprocess.Popen(['gmx', 'energy', '-f', edr, '-b', str(begin)], stdout=PIPE, stdin=PIPE, stderr=PIPE)
        sp_out = sp.communicate(input=property_str.encode())[0]

        for line in sp_out.decode().splitlines():
            if line.lower().startswith(property_str.lower()):
                return float(line.split()[1])
        raise Exception('Invalid property')

    @staticmethod
    def get_box(edr, begin=0) -> [float]:
        sp = subprocess.Popen(['gmx', 'energy', '-f', edr, '-b', str(begin)], stdout=PIPE, stdin=PIPE, stderr=PIPE)
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
    def scale_box(gro, new_gro, new_box:[float]):
        n = 0
        NAtoms = 0
        nresidue = []
        residue = []
        element = []
        natom = []
        xyz = []
        vxyz = []
        box = []
        if not os.path.exists(gro):
            gro = gro + '.gro'
        if not os.path.exists(gro):
            raise Exception('gro not found')
        with open(gro) as f:
            for line in f:
                n += 1
                if n == 1:
                    continue
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
                    xyz.append([x,y,z])
                    vxyz.append([vx,vy,vz])
                if n == 3 + NAtoms:
                    box = [ float(word) for word in line.strip().split()[:3]]
                    break

        scale = [new_box[i] / box[i] for i in range(3)]
        xyz = [[i[0] * scale[0], i[1]*scale[1], i[2]*scale[2]] for i in xyz]
        vxyz = [[i[0] * scale[0], i[1]*scale[1], i[2]*scale[2]] for i in vxyz]

        with open(new_gro, 'w') as f:
            f.write('Scaled box\n%i\n' %NAtoms)
            for i in range(NAtoms):
                f.write('%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n'
                        %(nresidue[i], residue[i], element[i], natom[i],
                          xyz[i][0],xyz[i][1],xyz[i][2],vxyz[i][0],vxyz[i][1],vxyz[i][2]))
            f.write('%f %f %f\n' %(new_box[0], new_box[1], new_box[2]))

    @staticmethod
    def generate_top(itp, molecules, numbers):
        shutil.copy(itp, 'topol.top')
        with open('topol.top', 'a') as f:
            f.write('\n[system]\n[molecules]\n')
            for i, molecule in enumerate(molecules):
                f.write('%s %i\n' %(molecule, numbers[i]))

    @staticmethod
    def pdb2gro(pdb, gro, box: [float]):
        if len(box) != 3:
            raise Exception('Invalid box')
        sp = subprocess.run(['gmx', 'editconf', '-f', pdb, '-o', gro, '-box', str(box[0]), str(box[1]), str(box[2])])


    @staticmethod
    def velacc(trr, tpr=None, group=None, begin=0, out='velacc', silent=False):
        if tpr == None:
            tpr = trr
        if group == None:
            raise Exception('No group specifed')

        if silent:
            sp = subprocess.Popen(['gmx', 'velacc', '-f', trr, '-s', tpr, '-o', out, '-b', str(begin), '-mol', '-nonormalize'],
                    stdin=PIPE, stdout=PIPE, stderr=PIPE)
        else:
            sp = subprocess.Popen(['gmx', 'velacc', '-f', trr, '-s', tpr, '-o', out, '-b', str(begin), '-mol', '-nonormalize'],
                    stdin=PIPE)
        sp.communicate(input=str(group).encode())

    @staticmethod
    def msd_com(xtc, tpr, resname, beginfit=-1, endfit=-1, out=None, silent=False):
        ndx = 'com-' + resname + '.ndx'
        GMX.select_com(tpr, resname, out=ndx)
        if out == None:
            out = 'msd-com-%s.xvg' %resname
        if silent:
            sp = subprocess.Popen(['gmx', 'msd', '-f', xtc, '-s', tpr, '-n', ndx, '-o', out, '-nomw', '-beginfit', str(beginfit), '-endfit', str(endfit)], stdout=PIPE, stderr=PIPE)
        else:
            sp = subprocess.Popen(['gmx', 'msd', '-f', xtc, '-s', tpr, '-n', ndx, '-o', out, '-nomw', '-beginfit', str(beginfit), '-endfit', str(endfit)])
        sp_out = sp.communicate()[0]

        for line in sp_out.decode().splitlines():
            if line.startswith('D['):
                return line
        raise Exception('Error running gmx msd')

    @staticmethod
    def traj_com(xtc, tpr, out=None, silent=False):
        ndx = 'com.ndx'
        GMX.select_com(tpr, 'all', out=ndx)
        if out == None:
            out = 'com.xtc'
        if silent:
            sp = subprocess.Popen(['gmx', 'traj', '-s', tpr, '-f', xtc, '-oxt', out, '-mol', '-n', ndx], stdout=PIPE, stderr=PIPE)
            sp.communicate()
        else:
            sp = subprocess.run(['gmx', 'traj', '-s', tpr, '-f', xtc, '-oxt', out, '-mol', '-n', ndx])
        return out, ndx

    @staticmethod
    def select_com(tpr, resname, out='selection.ndx'):
        sp = subprocess.Popen(['gmx', 'select', '-s', tpr, '-on', out], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        if resname == 'all':
            select_com_str = 'res_com of all'
        else:
            select_com_str = 'res_com of resname %s' %resname
        sp.communicate(input=select_com_str.encode())

