import math
from ..topology import Topology, Molecule, Atom


class Vector():
    """minimal 3-vector"""

    def __init__(self, x=0.0, y=0.0, z=0.0):
        if isinstance(x, tuple) or isinstance(x, list):
            self.x, self.y, self.z = x
        else:
            self.x = x
            self.y = y
            self.z = z

    def __getitem__(self, index):
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        else:
            raise IndexError('vector index out of range')

    def __abs__(self):
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def __add__(self, other):
        if isinstance(other, Vector):
            return Vector(self.x + other.x, self.y + other.y, self.z + other.z)
        else:
            raise TypeError('wrong type in vector addition')

    def __sub__(self, other):
        if isinstance(other, Vector):
            return Vector(self.x - other.x, self.y - other.y, self.z - other.z)
        else:
            raise TypeError('wrong type in vector subtraction')

    def __mul__(self, other):
        if isinstance(other, Vector):  # dot product
            return self.x * other.x + self.y * other.y + self.z * other.z
        else:
            return Vector(self.x * other, self.y * other, self.z * other)

    def __div__(self, other):
        return Vector(self.x / other, self.y / other, self.z / other)

    def __neg__(self):
        return Vector(-self.x, -self.y, -self.z)

    def __str__(self):
        return "( {0}, {1}, {2} )".format(self.x, self.y, self.z)

    def __repr__(self):
        return str(self) + ' instance at 0x' + str(hex(id(self))[2:].upper())

    def cross(self, other):
        return Vector(self.y * other.z - self.z * other.y,
                      self.z * other.x - self.x * other.z,
                      self.x * other.y - self.y * other.x)

    def unit(self):
        a = abs(self)
        return Vector(self.x / a, self.y / a, self.z / a)


class Zmat(Topology):
    """
    z-matrix representing a molecule, read from .zmat file
    TODO to be implemented
    """

    def __init__(self, file, mode='r'):
        super().__init__()
        self._file = open(file, mode)
        if mode == 'r':
            self.parse()
        elif mode == 'w':
            pass

    def __del__(self):
        self.close()

    def close(self):
        self._file.close()

    def parse(self):

        mol = Molecule()
        atoms = []
        zatom = []
        connect = []
        improper = []

        f = self._file

        # read molecule name
        line = f.readline()
        while line.strip().startswith('#'):
            line = f.readline()
        mol.name = line.strip()

        # read z-matrix
        line = f.readline()
        while line.strip().startswith('#') or line.strip() == '':
            line = f.readline()

        tok = line.strip().split()
        if len(tok) > 1:  # there can be line numbers
            shift = 1
        else:
            shift = 0

        variables = False
        while line and not line.strip().lower().startswith('var'):
            tok = line.strip().split()
            if len(tok) == 0:
                break
            name = tok[shift]
            ir = ia = id = 0
            r = a = d = 0.0
            rvar = avar = dvar = ''
            if (len(tok) - shift) > 1:
                ir = int(tok[shift + 1])
                if tok[shift + 2][0].isalpha():
                    rvar = tok[shift + 2]
                    variables = True
                else:
                    r = float(tok[shift + 2])
                if (len(tok) - shift) > 3:
                    ia = int(tok[shift + 3])
                    if tok[shift + 4][0].isalpha():
                        avar = tok[shift + 4]
                        variables = True
                    else:
                        a = float(tok[shift + 4])
                    if (len(tok) - shift) > 5:
                        id = int(tok[shift + 5])
                        if tok[shift + 6][0].isalpha():
                            dvar = tok[shift + 6]
                            variables = True
                        else:
                            d = float(tok[shift + 6])
            za = {'name': name,
                     'ir': ir, 'rvar': rvar, 'r': r,
                     'ia': ia, 'avar': avar, 'a': a,
                     'id': id, 'dvar': dvar, 'd': d}
            zatom.append(za)
            line = f.readline()

        # read variables
        if variables:
            if line.strip().lower().startswith('var') or line.strip() == '':
                line = f.readline()
            while line:
                tok = line.strip().split('=')
                if len(tok) < 2:
                    break
                key = tok[0].strip()
                val = float(tok[1])
                for rec in zatom:
                    if rec['rvar'] == key:
                        rec['r'] = val
                    if rec['avar'] == key:
                        rec['a'] = val
                    if rec['dvar'] == key:
                        rec['d'] = val
                line = f.readline()

        # read connects, improper, force field file
        self.ff = ''
        self.guessconnect = False
        while line:
            if line.strip().startswith('#') or line.strip() == '':
                line = f.readline()
                continue
            tok = line.strip().split()
            if tok[0] == 'reconnect':
                self.guessconnect = True
            if tok[0] == 'connect':
                atomi = int(tok[1])
                atomj = int(tok[2])
                connect.append([atomi, atomj])
            elif tok[0] == 'improper':
                atomi = int(tok[1])
                atomj = int(tok[2])
                atomk = int(tok[3])
                atoml = int(tok[4])
                improper.append([atomi, atomj, atomk, atoml])
            else:
                self.ff = tok[0]
            line = f.readline()

        natom = len(atoms)
        if natom != len(zatom):
            raise Exception('  error: different numbers of atoms in zmat ' + mol.name)

        if natom == 0:
            return self

        # first atom at origin
        atoms[0].x = 0.0
        atoms[0].y = 0.0
        atoms[0].z = 0.0
        if natom == 1:
            return self

        # second atom at distance r from first along xx
        atoms[1].x = zatom[1]['r']
        atoms[1].y = 0.0
        atoms[1].z = 0.0
        if natom == 2:
            return self

        # third atom at distance r from ir forms angle a 3-ir-ia in plane xy
        r = zatom[2]['r']
        ir = zatom[2]['ir'] - 1
        ang = zatom[2]['a'] * math.pi / 180.0
        ia = zatom[2]['ia'] - 1

        # for this construction, the new atom is at point (x, y), atom
        # ir is at point (xr, yr) and atom ia is at point (xa, ya).
        # Theta is the angle between the vector joining ir to ia and
        # the x-axis, a' (= theta - a) is is the angle between r and
        # the x-axis. x = xa + r cos a', y = ya + r sin a'.  From the
        # dot product of a unitary vector along x with the vector from
        # ir to ia, theta can be calculated: cos theta = (xa - xr) /
        # sqrt((xa - xr)^2 + (ya - yr)^2).  If atom ia is in third or
        # forth quadrant relative to atom ir, ya - yr < 0, then theta
        # = 2 pi - theta. */
        delx = atoms[ia].x - atoms[ir].x
        dely = atoms[ia].y - atoms[ir].y
        theta = math.acos(delx / math.sqrt(delx * delx + dely * dely))
        if dely < 0.0:
            theta = 2 * math.pi - theta
        ang = theta - ang
        atoms[2].x = atoms[ir].x + r * math.cos(ang)
        atoms[2].y = atoms[ir].y + r * math.sin(ang)
        atoms[2].z = 0.0
        if natom == 3:
            return self

        # nth atom at distance r from atom ir forms angle a at 3-ir-ia
        # and dihedral angle between planes 3-ir-ia and ir-ia-id
        for i in range(3, natom):
            r = zatom[i]['r']
            ir = zatom[i]['ir'] - 1
            ang = zatom[i]['a'] * math.pi / 180.0
            ia = zatom[i]['ia'] - 1
            dih = zatom[i]['d'] * math.pi / 180.0
            id = zatom[i]['id'] - 1

            # for this construction the new atom is at point A, atom ir is
            # at B, atom ia at C and atom id at D.  Point a is the
            # projection of A onto the plane BCD.  Point b is the
            # projection of A along the direction BC (the line defining
            # the dihedral angle between planes ABC and BCD). n = CD x BC
            # / |CD x BC| is the unit vector normal to the plane BCD. m =
            # BC x n / |BC x n| is the unit vector on the plane BCD normal
            # to the direction BC.
            #
            #                               .'A
            #                 ------------.' /.-----------------
            #                /           b /  .               /
            #               /           ./    .              /
            #              /           B......a      ^      /
            #             /           /              |n    /
            #            /           /                    /
            #           /           C                    /
            #          /             \                  /
            #         /               \                /
            #        /plane BCD        D              /
            #       ----------------------------------
            #
            #                    A              C------B...b
            #                   /.             /        .  .
            #                  / .            /    |m    . .
            #                 /  .           /     V      ..
            #         C------B...b          D              a
            #

            BA = r
            vB = Vector(atoms[ir].x, atoms[ir].y, atoms[ir].z)
            vC = Vector(atoms[ia].x, atoms[ia].y, atoms[ia].z)
            vD = Vector(atoms[id].x, atoms[id].y, atoms[id].z)

            vBC = vC - vB
            vCD = vD - vC

            BC = abs(vBC)
            bB = BA * math.cos(ang)
            bA = BA * math.sin(ang)
            aA = bA * math.sin(dih)
            ba = bA * math.cos(dih)

            vb = vC - vBC * ((BC - bB) / BC)
            vn = (vCD.cross(vBC)).unit()
            vm = (vBC.cross(vn)).unit()
            va = vb + vm * ba
            vA = va + vn * aA

            atoms[i].x = vA.x
            atoms[i].y = vA.y
            atoms[i].z = vA.z
