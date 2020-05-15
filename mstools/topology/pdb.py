import math
from ..topology import Topology
from ..forcefield import Element
from .atom import Atom
from .molecule import Molecule
from .topology import Topology


class Pdb(Topology):
    def __init__(self, file, **kwargs):
        super().__init__()
        self.parse(file)

    def parse(self, file):
        with open(file) as f:
            lines = f.read().splitlines()

        atoms = []
        mol_names = {}
        for line in lines:
            line = line.rstrip()
            if line.startswith('CRYST1'):
                a, b, c, alpha, beta, gamma = list(map(float, line.split()[1:7]))
                self.cell.set_box([[a / 10, b / 10, c / 10], [alpha, beta, gamma]])
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = Atom()
                atom.id = int(line[6:11]) - 1
                # only read the first frame
                if atom.id != len(atoms):
                    break
                atom.name = line[12:16].strip()
                atoms.append(atom)

                mol_name = line[17:21].strip()
                mol_id = int(line[22:26])
                if mol_id not in mol_names:
                    mol_names[mol_id] = mol_name
                atom._mol_id = mol_id  # temporary attribute for identifying molecules

                x = float(line[30:38]) / 10
                y = float(line[38:46]) / 10
                z = float(line[46:54]) / 10
                atom.position = (x, y, z,)
                if len(line) > 76:
                    symbol = line[76:78].strip()
                    element = Element(symbol)
                else:
                    element = Element.guess_from_atom_type(atom.name)
                atom.symbol = element.symbol
                atom.mass = element.mass

        molecules = [Molecule(name) for id, name in mol_names.items()]
        for atom in atoms:
            mol = molecules[atom._mol_id - 1]
            mol.add_atom(atom)

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
                atom1 = atoms[id1 - 1]
                atom2 = atoms[id2 - 1]
                atom1.molecule.add_bond(atom1, atom2, check_existence=True)
                if len(line) > 16:
                    id3 = int(line[16:21])
                    atom3 = atoms[id3 - 1]
                    atom1.molecule.add_bond(atom1, atom3, check_existence=True)
                if len(line) > 21:
                    id4 = int(line[21:26])
                    atom4 = atoms[id4 - 1]
                    atom1.molecule.add_bond(atom1, atom4, check_existence=True)
                if len(line) > 26:
                    id5 = int(line[26:31])
                    atom5 = atoms[id5 - 1]
                    atom1.molecule.add_bond(atom1, atom5, check_existence=True)
            prev_section = section

        self.update_molecules(molecules)

    @staticmethod
    def save_to(top, file):
        if not top.has_position:
            raise Exception('Position is required for writing PDB file')

        lengths = top.cell.lengths * 10
        angles = top.cell.angles
        with open(file, 'w') as f:
            f.write('REMARK   Created by mstools\n')
            f.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 \n' % (
                lengths[0], lengths[1], lengths[2], angles[0], angles[1], angles[2]))

            for atom in top.atoms:
                pos = atom.position * 10  # convert from nm to A
                line = 'HETATM%5d %4s %4s %4d    %8.3f%8.3f%8.3f                      %2s\n' % (
                    (atom.id + 1) % 100000, atom.name[:4], atom.molecule.name[:4],
                    (atom.molecule.id + 1) % 10000, pos[0], pos[1], pos[2], atom.symbol[:2])
                f.write(line)

            for atom in top.atoms:
                partners = [a for a in atom.bond_partners if a.id > atom.id]
                if len(partners) > 0:
                    line = 'CONECT%5d' % (atom.id + 1)
                    for i, partner in enumerate(partners):
                        if i > 0 and i % 4 == 0:
                            line += '\nCONECT%5d' % (atom.id + 1)
                        line += '%5d' % (partner.id + 1)
                    line += '\n'
                    f.write(line)
