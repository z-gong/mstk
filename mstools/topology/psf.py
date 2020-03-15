import math
from ..forcefield import Element
from ..topology import Atom, Molecule, Topology


class Psf(Topology):
    '''
    It parses psf file to generate Topology
    Atoms, bonds, angles, dihedrals and impropers are parsed,
    all others like hydrogen bonds are ignored
    Polarizability and thole parameters for Drude model are parsed,
    but anisotroic polarizability are ignored
    '''

    def __init__(self, file, mode='r'):
        super(Psf, self).__init__()
        self._file = open(file, mode)
        if mode == 'r':
            self.parse()
        elif mode == 'w':
            pass

    def parse(self):
        lines = self._file.read().splitlines()
        words = lines[0].strip().split()
        if 'PSF' not in words:
            raise Exception('Invalid psf file, PSF not found in first line')
        self.is_drude = 'DRUDE' in words

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
                self.parse_title(lines[iline: iline_next])
            if section == '!NATOM':
                self.parse_atoms(lines[iline: iline_next])
            if section == '!NBOND':
                self.parse_bonds(lines[iline: iline_next])
            if section == '!NTHETA':
                self.parse_angles(lines[iline: iline_next])
            if section == '!NPHI':
                self.parse_dihedrals(lines[iline: iline_next])
            if section == '!NIMPHI':
                self.parse_impropers(lines[iline: iline_next])

    def parse_title(self, lines):
        nline = int(lines[0].strip().split()[0])
        self.remark = '\n'.join(lines[1:1 + nline]).strip()

    def parse_atoms(self, lines):
        '''
        in psf file, drude particles always appear after attached atoms
        '''
        n_atom = int(lines[0].strip().split()[0])
        atoms = [Atom() for i in range(n_atom)]
        # i don't know how many molecules there will be
        # but it can not be larger than n_atom
        # i'll kick out the redundant molecules later
        molecules = [Molecule() for i in range(n_atom)]
        for i in range(n_atom):
            words = lines[1 + i].strip().split()
            atom_id = int(words[0])
            mol_id = int(words[2])
            mol_name = words[3]
            atom_type = words[5]
            charge, mass = list(map(float, words[6:8]))
            if self.is_drude:
                alpha, thole = list(map(float, words[9:11]))
            else:
                alpha = thole = 0.0

            mol = molecules[mol_id - 1]  # mol.id starts from 0
            if not mol.initialized:
                # set the information of this molecule
                mol.id = mol_id - 1
                mol.name = mol_name

            atom = atoms[atom_id - 1]  # atom.id starts from 0
            atom.id = atom_id - 1
            mol.add_atom(atom)
            atom.charge = charge
            atom.mass = mass
            atom.type = atom_type
            atom.symbol = Element.guess_from_atom_type(atom.type).symbol
            atom.alpha = -alpha
            atom.thole = thole * 2

        # kick out redundant molecules
        molecules = [mol for mol in molecules if mol.initialized]
        for mol in molecules:
            for i, atom in enumerate(mol.atoms):
                # atomic symbol + index inside mol starting from 1
                atom.name = atom.symbol + str(i + 1)

        self.init_from_molecules(molecules)

    def parse_bonds(self, lines):
        n_bond = int(lines[0].strip().split()[0])
        atoms = self.atoms
        for i in range(math.ceil(n_bond / 4)):
            words = lines[i + 1].strip().split()
            for j in range(4):
                if i * 4 + j + 1 > n_bond:
                    break
                atom1, atom2 = atoms[int(words[j * 2]) - 1], atoms[int(words[j * 2 + 1]) - 1]
                atom1.molecule.add_bond(atom1, atom2)

    def parse_angles(self, lines):
        n_angle = int(lines[0].strip().split()[0])
        atoms = self.atoms
        for i in range(math.ceil(n_angle / 3)):
            words = lines[i + 1].strip().split()
            for j in range(3):
                if i * 3 + j + 1 > n_angle:
                    break
                atom1, atom2, atom3 = atoms[int(words[j * 3 + 0]) - 1], \
                                      atoms[int(words[j * 3 + 1]) - 1], \
                                      atoms[int(words[j * 3 + 2]) - 1]
                atom1.molecule.add_angle(atom1, atom2, atom3)

    def parse_dihedrals(self, lines):
        n_dihedral = int(lines[0].strip().split()[0])
        atoms = self.atoms
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

    def parse_impropers(self, lines):
        n_improper = int(lines[0].strip().split()[0])
        atoms = self.atoms
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

    def write(self):
        if self.is_drude:
            self._file.write('PSF DRUDE\n\n')
        else:
            self._file.write('PSF\n\n')

        self._file.write('%8i !NTITLE\n' % 1)
        self._file.write(' REMARKS Created by mstools\n\n')

        self._file.write('%8i !NATOM\n' % self.n_atom)
        for i, atom in enumerate(self.atoms):
            line = '%8i S    %-4i %4s %-4s %-4s %10.6f %13.4f %11i %13.4f %13.4f\n' % (
                atom.id + 1, atom.molecule.id + 1, atom.molecule.name[:4], atom.symbol, atom.type,
                atom.charge, atom.mass, 0, -atom.alpha, atom.thole / 2
            )
            self._file.write(line)
        self._file.write('\n')

        n_bond = self.n_bond
        self._file.write('%8i !NBOND: bonds\n' % n_bond)
        for i, bond in enumerate(self.bonds):
            self._file.write('%8i%8i' % (bond.atom1.id + 1, bond.atom2.id + 1))
            if (i + 1) % 4 == 0 or i + 1 == n_bond:
                self._file.write('\n')
        self._file.write('\n')

        n_angle = self.n_angle
        self._file.write('%8i !NTHETA: angles\n' % n_angle)
        for i, angle in enumerate(self.angles):
            self._file.write('%8i%8i%8i' % (angle.atom1.id + 1, angle.atom2.id + 1, angle.atom3.id + 1))
            if (i + 1) % 3 == 0 or i + 1 == n_angle:
                self._file.write('\n')
        self._file.write('\n')

        n_dihedral = self.n_dihedral
        self._file.write('%8i !NPHI: dihedrals\n' % n_dihedral)
        for i, dihedral in enumerate(self.dihedrals):
            self._file.write('%8i%8i%8i%8i' % (dihedral.atom1.id + 1, dihedral.atom2.id + 1,
                                               dihedral.atom3.id + 1, dihedral.atom4.id + 1))
            if (i + 1) % 2 == 0 or i + 1 == n_dihedral:
                self._file.write('\n')
        self._file.write('\n')

        n_improper = self.n_improper
        self._file.write('%8i !NIMPHI: impropers\n' % n_improper)
        for i, improper in enumerate(self.impropers):
            self._file.write('%8i%8i%8i%8i' % (improper.atom1.id + 1, improper.atom2.id + 1,
                                               improper.atom3.id + 1, improper.atom4.id + 1))
            if (i + 1) % 2 == 0 or i + 1 == n_improper:
                self._file.write('\n')
        self._file.write('\n')

        self._file.write('%8i !NDON: donors\n\n' % 0)

        self._file.write('%8i !NACC: acceptors\n\n' % 0)

        self._file.write('%8i !NNB\n\n\n' % 0)

        self._file.write('%8i !NUMANISO\n\n' % 0)

        self._file.flush()
