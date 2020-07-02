import math
from ..forcefield import Element
from .atom import Atom
from .molecule import Molecule
from .topology import Topology


class Psf():
    '''
    Generate Topology from PSF file.

    Atoms, bonds, angles, dihedrals and impropers are parsed.
    Polarizability and thole parameters for Drude model are parsed,
    but anisotroic polarizability are ignored
    All others like hydrogen bonds are ignored

    This parser is designed for non-biological system.
    Each residue is treated as one molecule.
    So if there are bonds between residues, there might be bugs.

    * Todo Support inter-residue bonds

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
    >>> psf = Psf('input.psf')
    >>> topology = psf.topology

    >>> Psf.save_to(topology, 'output.psf')
    '''

    def __init__(self, file, **kwargs):
        self.topology = Topology()
        self._parse(file)

    def _parse(self, file):
        with open(file) as f:
            lines = f.read().splitlines()

        words = lines[0].strip().split()
        if 'PSF' not in words:
            raise Exception('Invalid psf file, PSF not found in first line')
        is_drude = 'DRUDE' in words

        # there are other sections, but i don't care
        _sections = ['!NTITLE', '!NATOM', '!NBOND', '!NTHETA', '!NPHI', '!NIMPHI']

        # record the line number (starts from 0) of section headers
        _iline_section = []
        for i, line in enumerate(lines):
            words = line.strip().split()
            if len(words) < 2:
                continue
            word = words[1].split(':')[0]
            if word in _sections:
                _iline_section.append((i, word))

        for i, (iline, section) in enumerate(_iline_section):
            if i != len(_iline_section) - 1:
                iline_next = _iline_section[i + 1][0]
            else:
                iline_next = len(lines)
            if section == '!NTITLE':
                # the number at section line indicates the line counts of the section
                self._parse_title(lines[iline: iline_next])
            if section == '!NATOM':
                self._parse_atoms(lines[iline: iline_next], is_drude)
            if section == '!NBOND':
                self._parse_bonds( lines[iline: iline_next])
            if section == '!NTHETA':
                self._parse_angles(lines[iline: iline_next])
            if section == '!NPHI':
                self._parse_dihedrals(lines[iline: iline_next])
            if section == '!NIMPHI':
                self._parse_impropers(lines[iline: iline_next])

    def _parse_title(self, lines):
        nline = int(lines[0].strip().split()[0])
        self.topology.remark = '\n'.join(lines[1:1 + nline]).strip()

    def _parse_atoms(self, lines, parse_drude):
        '''
        in psf file, drude particles always appear after attached atoms
        '''
        n_atom = int(lines[0].strip().split()[0])
        atoms = [Atom() for i in range(n_atom)]
        mol_names = {}
        for i in range(n_atom):
            words = lines[1 + i].strip().split()
            atom_id = int(words[0])
            mol_id = int(words[2])
            mol_name = words[3]
            atom_name = words[4]
            atom_type = words[5]
            charge, mass = list(map(float, words[6:8]))
            if parse_drude:
                alpha, thole = list(map(float, words[9:11]))
            else:
                alpha = thole = 0.0

            if mol_id not in mol_names:
                mol_names[mol_id] = mol_name

            atom = atoms[atom_id - 1]  # atom.id starts from 0
            atom.id = atom_id - 1
            atom.name = atom_name
            atom.charge = charge
            atom.mass = mass
            atom.type = atom_type
            atom._mol_id = mol_id  # temporary attribute for identifying molecules
            atom.alpha = -alpha / 1000  # convert A^3 to nm^3
            atom.thole = thole * 2
            if parse_drude and atom.type.startswith('D'):
                atom.is_drude = True
                atom.symbol = 'DP'
            else:
                atom.symbol = Element.guess_from_atom_type(atom.type).symbol

        molecules = [Molecule(name) for id, name in mol_names.items()]
        for atom in atoms:
            mol = molecules[atom._mol_id - 1]
            mol.add_atom(atom)

        self.topology.update_molecules(molecules)

    def _parse_bonds(self, lines):
        n_bond = int(lines[0].strip().split()[0])
        atoms = self.topology.atoms
        for i in range(math.ceil(n_bond / 4)):
            words = lines[i + 1].strip().split()
            for j in range(4):
                if i * 4 + j + 1 > n_bond:
                    break
                atom1, atom2 = atoms[int(words[j * 2]) - 1], atoms[int(words[j * 2 + 1]) - 1]
                atom1.molecule.add_bond(atom1, atom2)

    def _parse_angles(self, lines):
        n_angle = int(lines[0].strip().split()[0])
        atoms = self.topology.atoms
        for i in range(math.ceil(n_angle / 3)):
            words = lines[i + 1].strip().split()
            for j in range(3):
                if i * 3 + j + 1 > n_angle:
                    break
                atom1, atom2, atom3 = atoms[int(words[j * 3 + 0]) - 1], \
                                      atoms[int(words[j * 3 + 1]) - 1], \
                                      atoms[int(words[j * 3 + 2]) - 1]
                atom1.molecule.add_angle(atom1, atom2, atom3)

    def _parse_dihedrals(self, lines):
        n_dihedral = int(lines[0].strip().split()[0])
        atoms = self.topology.atoms
        for i in range(math.ceil(n_dihedral / 2)):
            words = lines[i + 1].strip().split()
            for j in range(2):
                if i * 2 + j + 1 > n_dihedral:
                    break
                atom1, atom2, atom3, atom4 = atoms[int(words[j * 4 + 0]) - 1], \
                                             atoms[int(words[j * 4 + 1]) - 1], \
                                             atoms[int(words[j * 4 + 2]) - 1], \
                                             atoms[int(words[j * 4 + 3]) - 1]
                atom1.molecule.add_dihedral(atom1, atom2, atom3, atom4)

    def _parse_impropers(self, lines):
        n_improper = int(lines[0].strip().split()[0])
        atoms = self.topology.atoms
        for i in range(math.ceil(n_improper / 2)):
            words = lines[i + 1].strip().split()
            for j in range(2):
                if i * 2 + j + 1 > n_improper:
                    break
                atom1, atom2, atom3, atom4 = atoms[int(words[j * 4 + 0]) - 1], \
                                             atoms[int(words[j * 4 + 1]) - 1], \
                                             atoms[int(words[j * 4 + 2]) - 1], \
                                             atoms[int(words[j * 4 + 3]) - 1]
                atom1.molecule.add_improper(atom1, atom2, atom3, atom4)

    @staticmethod
    def save_to(top, file):
        '''
        Save topology into a PSF file

        Parameters
        ----------
        top : Topology
        file : str
        '''
        for atom in top.atoms:
            if atom.type == '':
                raise Exception('PSF requires all atom types are defined')

        f = open(file, 'w')

        if top.is_drude:
            f.write('PSF DRUDE\n\n')
        else:
            f.write('PSF\n\n')

        f.write('%8i !NTITLE\n' % 1)
        f.write(' REMARKS Created by mstools\n\n')

        f.write('%8i !NATOM\n' % top.n_atom)
        for i, atom in enumerate(top.atoms):
            line = '%7i  S  %4i %8s %8s %8s %10.6f %8.4f %4i %8.4f %8.4f\n' % (
                atom.id + 1, atom.molecule.id + 1, atom.molecule.name, atom.name, atom.type,
                atom.charge, atom.mass, 0, -atom.alpha * 1000, atom.thole / 2
            )
            f.write(line)
        f.write('\n')

        n_bond = top.n_bond
        f.write('%8i !NBOND: bonds\n' % n_bond)
        for i, bond in enumerate(top.bonds):
            f.write('%8i%8i' % (bond.atom1.id + 1, bond.atom2.id + 1))
            if (i + 1) % 4 == 0 or i + 1 == n_bond:
                f.write('\n')
        f.write('\n')

        n_angle = top.n_angle
        f.write('%8i !NTHETA: angles\n' % n_angle)
        for i, angle in enumerate(top.angles):
            f.write('%8i%8i%8i' % (angle.atom1.id + 1, angle.atom2.id + 1,
                                   angle.atom3.id + 1))
            if (i + 1) % 3 == 0 or i + 1 == n_angle:
                f.write('\n')
        f.write('\n')

        n_dihedral = top.n_dihedral
        f.write('%8i !NPHI: dihedrals\n' % n_dihedral)
        for i, dihedral in enumerate(top.dihedrals):
            f.write('%8i%8i%8i%8i' % (dihedral.atom1.id + 1, dihedral.atom2.id + 1,
                                      dihedral.atom3.id + 1, dihedral.atom4.id + 1))
            if (i + 1) % 2 == 0 or i + 1 == n_dihedral:
                f.write('\n')
        f.write('\n')

        n_improper = top.n_improper
        f.write('%8i !NIMPHI: impropers\n' % n_improper)
        for i, improper in enumerate(top.impropers):
            f.write('%8i%8i%8i%8i' % (improper.atom1.id + 1, improper.atom2.id + 1,
                                      improper.atom3.id + 1, improper.atom4.id + 1))
            if (i + 1) % 2 == 0 or i + 1 == n_improper:
                f.write('\n')
        f.write('\n')

        f.write('%8i !NDON: donors\n\n' % 0)
        f.write('%8i !NACC: acceptors\n\n' % 0)
        f.write('%8i !NNB\n\n\n' % 0)
        f.write('%8i !NUMANISO\n\n' % 0)
        f.close()
