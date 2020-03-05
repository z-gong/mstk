from . import Element, Atom, Molecule, Topology

class LammpsData(Topology):
    '''
    Parse the topology information from Lammps data file
    Currently only atom type, charge and molecule are parsed
    The positions and connectivity haven't been implemented
    '''
    def __init__(self, file, mode='r'):
        super(LammpsData, self).__init__()
        self._file = open(file)
        if mode == 'r':
            self.parse()
        if mode == 'w':
            raise Exception('Writing support for LammpsData haven\'t been implemented')

    def parse(self):
        lines = self._file.read().splitlines()
        self.remark = lines[0].strip()

        _sections = ["Atoms", "Velocities", "Ellipsoids", "Lines", "Triangles", "Bodies",
                     "Bonds", "Angles", "Dihedrals", "Impropers",
                     "Masses", "Pair Coeffs", "PairIJ Coeffs", "Bond Coeffs", "Angle Coeffs",
                     "Dihedral Coeffs", "Improper Coeffs",
                     "BondBond Coeffs", "BondAngle Coeffs", "MiddleBondTorsion Coeffs",
                     "EndBondTorsion Coeffs", "AngleTorsion Coeffs",
                     "AngleAngleTorsion Coeffs", "BondBond13 Coeffs", "AngleAngle Coeffs"]

        # record the line number (starts from 0) of section headers
        _iline_section = []
        for i, line in enumerate(lines):
            line_clean = line.split('#')[0].strip()
            if line_clean in _sections:
                _iline_section.append((i, line_clean))

        self.parse_header(lines[:_iline_section[0][0]])

        self.atoms = [Atom() for i in range(self.n_atom)]
        self._type_masses = [0.] * (self.n_atom_type + 1)
        self._type_names = [''] * (self.n_atom_type + 1)
        self._molecule_dict = {}  # because the length of molecules is unknown. So use dict instead of list

        for i, (iline, section) in enumerate(_iline_section):
            if i != len(_iline_section) - 1:
                iline_next = _iline_section[i + 1][0]
            else:
                iline_next = len(lines)
            if section == 'Masses':
                self.parse_masses(lines[iline + 2: iline_next])  # there is a blank line after each section header
            if section == 'Atoms':
                self.parse_atoms(lines[iline + 2: iline_next])
            if section == 'Bonds':
                self.parse_bonds(lines[iline + 2: iline_next])

    def parse_header(self, lines):
        for line in lines:
            if line.endswith(' atoms'):
                self.n_atom = int(line.split()[0])
            elif line.endswith(' bonds'):
                self.n_bond = int(line.split()[0])
            elif line.endswith(' angles'):
                self.n_angle = int(line.split()[0])
            elif line.endswith(' dihedrals'):
                self.n_dihedral = int(line.split()[0])
            elif line.endswith(' impropers'):
                self.n_improper = int(line.split()[0])
            elif line.endswith(' atom types'):
                self.n_atom_type = int(line.split()[0])
            elif line.endswith(' bond types'):
                self.n_bond_type = int(line.split()[0])
            elif line.endswith(' angle types'):
                self.n_angle_type = int(line.split()[0])
            elif line.endswith(' dihedral types'):
                self.n_dihedral_type = int(line.split()[0])
            elif line.endswith(' improper types'):
                self.n_improper_type = int(line.split()[0])
            else:
                continue

    def parse_masses(self, lines):
        for i in range(self.n_atom_type):
            words = lines[i].strip().split()
            type_id = int(words[0])
            self._type_masses[type_id] = float(words[1])
            if len(words) > 2 and words[2] == '#':
                if words[-1] == 'DP':
                    self._type_names[type_id] = 'D'  # Drude particles
                else:
                    self._type_names[type_id] = words[3]
            else:
                self._type_names[type_id] = 'AT' + str(type_id)

    def parse_atoms(self, lines):
        for i in range(self.n_atom):
            words = lines[i].strip().split()
            atom_id = int(words[0])
            mol_id = int(words[1])
            type_id = int(words[2])
            charge = float(words[3])

            self.n_molecule = max(self.n_molecule, mol_id)
            if mol_id not in self._molecule_dict.keys():
                mol_name = words[9] if (len(words) > 7 and words[7] == '#') else 'UNK'
                mol = Molecule(mol_name)
                mol.id = mol_id - 1  # mol.id starts from 0
                self._molecule_dict[mol_id] = mol
            else:
                mol = self._molecule_dict[mol_id]

            atom = self.atoms[atom_id - 1]
            atom.id = atom_id - 1  # atom.id starts from 0
            mol.add_atom(atom)
            atom.charge = charge
            atom.mass = self._type_masses[type_id]
            atom.type = self._type_names[type_id]
            atom.symbol = Element.guess_from_atom_type(atom.type).symbol

        self.molecules = [self._molecule_dict[i + 1] for i in range(self.n_molecule)]
        for mol in self.molecules:
            for i, atom in enumerate(mol.atoms):
                atom.name = atom.symbol + str(i + 1)  # atomic symbol + index inside mol starting from 1

    def parse_bonds(self, lines):
        pass

    def parse_angles(self, lines):
        pass

    def parse_dihedrals(self, lines):
        pass

    def parse_impropers(self, lines):
        pass


