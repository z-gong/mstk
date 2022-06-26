from mstk.chem.element import Element
from mstk.topology.atom import Atom
from mstk.topology.molecule import Molecule
from mstk.topology.topology import Topology


class Pdb():
    '''
    Generate Topology from PDB file.

    Atom symbols, molecule information, unit cell and positions are parsed.
    The connectivity is parsed if CONECT section is provided in the PDB file.

    This parser is designed for non-biological system.
    Each residue is treated as one molecule.
    So if there are bonds between residues, there might be bugs.

    Parameters
    ----------
    file : str
    atom_type : bool
        If True, the atom name column will be parsed as atom type

    Attributes
    ----------
    topology : Topology

    Examples
    --------
    >>> pdb = Pdb('input.pdb')
    >>> topology = pdb.topology

    >>> Pdb.save_to(topology, 'output.pdb')
    '''

    def __init__(self, file, **kwargs):
        self.topology = Topology()
        self._molecule = Molecule()
        self._parse(file, **kwargs)

    def _parse(self, file, atom_type=False):
        '''
        Parse a PDB file


        Returns
        -------
        topology : Topology
        '''

        with open(file) as f:
            lines = f.read().splitlines()

        mol = self._molecule
        atom_ids = {}  # {int: Atom}
        res_names = {}  # {int: str}
        res_atoms = {}  # {int: [Atom]}
        _atom_started = False
        for line in lines:
            line = line.rstrip()
            if _atom_started and line == 'ENDMDL':
                break
            if line.startswith('CRYST1'):
                a, b, c, alpha, beta, gamma = list(map(float, line.split()[1:7]))
                self.topology.cell.set_box([[a / 10, b / 10, c / 10], [alpha, beta, gamma]])
            if line.startswith('ATOM') or line.startswith('HETATM'):
                _atom_started = True
                atom_id = int(line[6:11])
                atom_name = line[12:16].strip()
                res_name = line[17:21].strip()
                res_id = int(line[22:26])
                x = float(line[30:38]) / 10
                y = float(line[38:46]) / 10
                z = float(line[46:54]) / 10
                if len(line) > 76:
                    symbol = line[76:78].strip()
                    element = Element(symbol)
                else:
                    element = Element.guess_from_atom_type(atom_name)

                atom = Atom(atom_name)
                atom.position = (x, y, z,)
                atom.symbol = element.symbol
                atom.mass = element.mass
                if atom_type:
                    atom.type = atom_name
                else:
                    atom.name = atom_name

                mol.add_atom(atom)

                atom_ids[atom_id] = atom
                if res_id not in res_names:
                    res_names[res_id] = res_name
                    res_atoms[res_id] = []
                res_atoms[res_id].append(atom)

        for resid, resname in res_names.items():
            atoms = res_atoms[resid]
            mol.add_residue(resname, atoms)

        prev_section = ''
        for line in lines:
            line = line.rstrip()
            section = line[:6]
            # only read the first frame
            if prev_section == 'CONECT' and section != 'CONECT':
                break
            if section == 'CONECT':
                id1 = int(line[6:11])
                id2 = int(line[11:16])
                atom1 = atom_ids[id1]
                atom2 = atom_ids[id2]
                mol.add_bond(atom1, atom2, check_existence=True)
                if len(line) > 16:
                    id3 = int(line[16:21])
                    atom3 = atom_ids[id3]
                    mol.add_bond(atom1, atom3, check_existence=True)
                if len(line) > 21:
                    id4 = int(line[21:26])
                    atom4 = atom_ids[id4]
                    mol.add_bond(atom1, atom4, check_existence=True)
                if len(line) > 26:
                    id5 = int(line[26:31])
                    atom5 = atom_ids[id5]
                    mol.add_bond(atom1, atom5, check_existence=True)
            prev_section = section

        self.topology.update_molecules(mol.split(consecutive=True))
        self.topology.generate_angle_dihedral_improper()

    @staticmethod
    def save_to(top, file, atom_type=False, **kwargs):
        '''
        Save topology into a PDB file

        Parameters
        ----------
        top : Topology
        file : str
        atom_type : bool
            If True, atom type will be written in the atom name column, which is useful for visualizing
        '''
        if not top.has_position:
            raise Exception('Position is required for writing PDB file')

        lengths = top.cell.lengths * 10
        angles = top.cell.angles
        string = 'REMARK   Created by mstk\n'
        string += 'CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 \n' % (
            lengths[0], lengths[1], lengths[2], angles[0], angles[1], angles[2])

        for atom in top.atoms:
            pos = atom.position * 10  # convert from nm to A
            atom_name = atom.type if atom_type else atom.name
            resname = atom.residue.name
            resid = atom.residue.id + 1
            line = 'HETATM%5d %4s %4s %4d    %8.3f%8.3f%8.3f                      %2s\n' % (
                (atom.id + 1) % 100000, atom_name[:4], resname[:4], resid % 10000,
                pos[0], pos[1], pos[2], atom.symbol[:2])
            string += line

        for atom in top.atoms:
            partners = [a for a in atom.bond_partners if a.id > atom.id]
            if len(partners) > 0:
                line = 'CONECT%5d' % (atom.id + 1)
                for i, partner in enumerate(partners):
                    if i > 0 and i % 4 == 0:
                        line += '\nCONECT%5d' % (atom.id + 1)
                    line += '%5d' % (partner.id + 1)
                line += '\n'
                string += line

        with open(file, 'wb') as f:
            f.write(string.encode())


Topology.register_format('.pdb', Pdb)
