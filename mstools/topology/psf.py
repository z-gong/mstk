from ..topology import Atom, Molecule, Topology


class Psf(Topology):
    def __init__(self, file, mode='r'):
        super(Psf, self).__init__()
        self._file = open(file, mode)
        if mode == 'r':
            self.parse()
        elif mode == 'w':
            raise Exception('Writing support for PSF haven\'t been implemented')

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
                self.parse_title(lines[iline: iline_next])  # the number at section line indicates the line counts of the section
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
            words = lines[1+i].strip().split()
            atom_id = int(words[0])
            mol_id = int(words[2])
            mol_name = words[3]
            atom_symbol = words[4]
            atom_type = words[5]
            charge, mass = list(map(float, words[6:8]))
            if self.is_drude:
                # ignore for now
                alpha, thole = list(map(float, words[9:11]))

            mol = molecules[mol_id - 1] # mol.id starts from 0
            if not mol.initiated:
                # set the information of this molecule
                mol.id = mol_id - 1
                mol.name = mol_name

            atom = atoms[atom_id - 1] # atom.id starts from 0
            atom.id = atom_id - 1
            mol.add_atom(atom)
            atom.charge = charge
            atom.mass = mass
            atom.type = atom_type
            atom.symbol = atom_symbol

        # kick out redundant molecules
        molecules = [mol for mol in molecules if mol.initiated]
        for mol in molecules:
            for i, atom in enumerate(mol.atoms):
                atom.name = atom.symbol + str(i + 1)  # atomic symbol + index inside mol starting from 1

        self.init_from_molecules(molecules)

    def parse_bonds(self, lines):
        pass

    def parse_angles(self, lines):
        pass

    def parse_dihedrals(self, lines):
        pass

    def parse_impropers(self, lines):
        pass
