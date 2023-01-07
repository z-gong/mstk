import math
import numpy as np
from mstk.topology.atom import Atom
from mstk.topology.molecule import Molecule
from mstk.topology.topology import Topology
from mstk.chem.element import Element
from mstk.chem.constant import *


class Zmat():
    '''
    Generate Topology from ZMAT file.

    All of the atoms are assumed to be in the same molecule.
    The first line is treated as the name of the molecule.
    The positions will be generated from Z-matrix, with the algorithm copied from fftool of Agilio Padua.
    The first column is treated as atom type.

    Parameters
    ----------
    file : str
    kwargs : dict
        Ignored

    Attributes
    ----------
    topology : Topology

    Examples
    --------
    >>> zmat = Zmat('input.zmat')
    >>> topology = zmat.topology
    '''

    def __init__(self, file, **kwargs):
        self.topology = Topology()
        self._parse(file)

    def _parse(self, file):
        f = open(file)

        mol = Molecule()
        z_atoms = []

        # read molecule name
        line = f.readline()
        while line.strip().startswith('#'):
            line = f.readline()
        mol.name = line.strip()

        # read z-matrix
        line = f.readline()
        while line.strip().startswith('#') or line.strip() == '':
            line = f.readline()

        words = line.strip().split()
        if len(words) > 1:  # there can be line numbers
            shift = 1
        else:
            shift = 0

        _read_variable = False
        while line and not line.strip().lower().startswith('var'):
            words = line.strip().split()
            if len(words) == 0:
                break
            atom = Atom()
            atom.type = words[shift]
            element = Element.guess_from_atom_type(atom.type)
            atom.symbol = element.symbol
            atom.mass = element.mass
            atom.name = atom.symbol + str(mol.n_atom + 1)
            mol.add_atom(atom)

            ir = ia = id = 0
            r = a = d = 0.0
            var_r = var_a = var_d = ''
            if (len(words) - shift) > 1:
                ir = int(words[shift + 1]) - 1
                mol.add_bond(mol.atoms[ir], atom)

                if words[shift + 2][0].isalpha():
                    var_r = words[shift + 2]
                    _read_variable = True
                else:
                    r = float(words[shift + 2]) / 10  # convert from A to nm
                if (len(words) - shift) > 3:
                    ia = int(words[shift + 3]) - 1
                    if words[shift + 4][0].isalpha():
                        var_a = words[shift + 4]
                        _read_variable = True
                    else:
                        a = float(words[shift + 4])
                    if (len(words) - shift) > 5:
                        id = int(words[shift + 5]) - 1
                        if words[shift + 6][0].isalpha():
                            var_d = words[shift + 6]
                            _read_variable = True
                        else:
                            d = float(words[shift + 6])
            za = {'name': type,
                  'ir'  : ir, 'var_r': var_r, 'r': r,
                  'ia'  : ia, 'var_a': var_a, 'a': a,
                  'id'  : id, 'var_d': var_d, 'd': d}
            z_atoms.append(za)
            line = f.readline()

        # read variables
        if _read_variable:
            if line.strip().lower().startswith('var') or line.strip() == '':
                line = f.readline()
            while line:
                words = line.strip().split('=')
                if len(words) < 2:
                    break
                key = words[0].strip()
                val = float(words[1])
                for za in z_atoms:
                    if za['var_r'] == key:
                        za['r'] = val / 10  # covert from A to nm
                    if za['var_a'] == key:
                        za['a'] = val
                    if za['var_d'] == key:
                        za['d'] = val
                line = f.readline()

        # read extra bonds
        while line:
            if line.strip().startswith('#') or line.strip() == '':
                line = f.readline()
                continue
            words = line.strip().split()
            if words[0] == 'connect':
                ii = int(words[1]) - 1
                ij = int(words[2]) - 1
                mol.add_bond(mol.atoms[ii], mol.atoms[ij])
            line = f.readline()

        f.close()

        if mol.n_atom >= 1:
            mol.atoms[0].position = [0, 0, 0]

        if mol.n_atom >= 2:
            mol.atoms[1].position = [z_atoms[1]['r'], 0, 0]

        # third atom at distance r from ir forms angle a 3-ir-ia in plane xy
        if mol.n_atom >= 3:
            r = z_atoms[2]['r']
            ir = z_atoms[2]['ir']
            ang = z_atoms[2]['a'] * DEG2RAD
            ia = z_atoms[2]['ia']

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
            dx, dy, dz = mol.atoms[ia].position - mol.atoms[ir].position
            theta = math.acos(np.clip(dx / math.sqrt(dx * dx + dy * dy), -1, 1))
            if dy < 0.0:
                theta = 2 * math.pi - theta
            ang = theta - ang
            mol.atoms[2].position = [mol.atoms[ir].position[0] + r * math.cos(ang),
                                     mol.atoms[ir].position[1] + r * math.sin(ang),
                                     0.0]

        # nth atom at distance r from atom ir forms angle a at 3-ir-ia
        # and dihedral angle between planes 3-ir-ia and ir-ia-id
        if mol.n_atom >= 4:
            for i in range(3, mol.n_atom):
                r = z_atoms[i]['r']
                ir = z_atoms[i]['ir']
                ang = z_atoms[i]['a'] * DEG2RAD
                ia = z_atoms[i]['ia']
                dih = z_atoms[i]['d'] * DEG2RAD
                id = z_atoms[i]['id']

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
                vB = np.array(mol.atoms[ir].position)
                vC = np.array(mol.atoms[ia].position)
                vD = np.array(mol.atoms[id].position)

                vBC = vC - vB
                vCD = vD - vC

                BC = math.sqrt(vBC.dot(vBC))
                bB = BA * math.cos(ang)
                bA = BA * math.sin(ang)
                aA = bA * math.sin(dih)
                ba = bA * math.cos(dih)

                vb = vC - vBC * ((BC - bB) / BC)
                vn = np.cross(vCD, vBC)
                vn = vn / math.sqrt(vn.dot(vn))
                vm = np.cross(vBC, vn)
                vm = vm / math.sqrt(vm.dot(vm))
                va = vb + vm * ba
                vA = va + vn * aA

                mol.atoms[i].position = vA

        mol.generate_angle_dihedral_improper()

        self.topology.update_molecules([mol])


Topology.register_format('.zmat', Zmat)
