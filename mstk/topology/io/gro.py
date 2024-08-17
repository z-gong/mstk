import numpy as np
from mstk.chem.element import Element
from mstk.topology.atom import Atom
from mstk.topology.molecule import Molecule
from mstk.topology.topology import Topology


class GroTopology:
    '''
    Generate Topology from GRO file.

    Atom name, residue information, unit cell and positions are parsed.
    The atom name colume will be used as both atom name and atom type.
    GRO file does not contain connectivity information. Each residue will be considered as a molecule.

    Only the first frame will be considered.

    Parameters
    ----------
    file : str

    Attributes
    ----------
    topology : Topology

    Examples
    --------
    >>> gro = Gro('input.gro')
    >>> topology = gro.topology

    >>> Gro.save_to(topology, 'output.gro')
    '''

    def __init__(self, file, **kwargs):
        self.topology = Topology()
        self._parse(file)

    def _parse(self, file):
        molecules = {}  # {int: Molecule}
        n_atom = 0
        last_resid = None
        with open(file) as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                if i == 1:
                    n_atom = int(line)
                    continue
                if i == n_atom + 2:
                    _box = tuple(map(float, line.split()))
                    if len(_box) == 3:
                        self.topology.cell.set_box(_box)
                    elif len(_box) == 9:
                        ax, by, cz, ay, az, bx, bz, cx, cy = _box
                        self.topology.cell.set_box([[ax, ay, az], [bx, by, bz], [cx, cy, cz]])
                    else:
                        raise ValueError('Invalid box')
                    break

                res_id = int(line[:5])
                res_name = line[5:10].strip()
                atom_name = line[10:15].strip()
                element = Element.guess_from_atom_type(atom_name)
                atom_id = int(line[15:20])
                x = float(line[20:28])
                y = float(line[28:36])
                z = float(line[36:44])

                atom = Atom(atom_name)
                atom.type = atom_name
                atom.symbol = element.symbol
                atom.mass = element.mass
                atom.position = (x, y, z,)

                if res_id != last_resid:
                    mol = Molecule(res_name)
                    molecules[res_id] = mol
                    last_resid = res_id
                else:
                    mol = molecules[res_id]
                mol.add_atom(atom)

        self.topology.update_molecules(list(molecules.values()))

    @staticmethod
    def save_to(top, file, **kwargs):
        '''
        Save topology into a GRO file

        Parameters
        ----------
        top : Topology
        file : str
        '''
        if not top.has_position:
            raise Exception('Position is required for writing GRO file')

        def format_float(value, max_length, default_decimal):
            if value >= 10 ** max_length or value <= -10 ** (max_length - 1):
                raise Exception(f'{value} too long for formatter')
            formatter = f'%{max_length}.{default_decimal}f'
            return (formatter % value)[:max_length]

        string = 'Created by mstk\n'
        string += f'{top.n_atom}\n'

        for atom in top.atoms:
            pos = atom.position
            if any(np.abs(pos) >= 1000):
                raise Exception('Positions are too large to be written in GRO format')
            string += "%5i%-5s%5s%5i%8.3f%8.3f%8.3f\n" % (
                (atom.residue.id + 1) % 100000, atom.residue.name[:5], atom.name[:5], (atom.id + 1) % 100000,
                pos[0], pos[1], pos[2]
            )

        if top.cell.is_rectangular:
            string += ' %.4f %.4f %.4f\n' % tuple(top.cell.get_size())
        else:
            a, b, c = top.cell.vectors
            string += ' %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n' % (
                a[0], b[1], c[2], a[1], a[2], b[0], b[2], c[0], c[1])

        with open(file, 'wb') as f:
            f.write(string.encode())


Topology.register_format('.gro', GroTopology)
