from ..forcefield import Element
from ..topology import Atom, Molecule, Topology

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

        n_atom, n_bond, n_angle, n_dihedral, n_improper, n_atom_type = self.parse_header(lines[:_iline_section[0][0]])

        self._type_masses = [0.] * (n_atom_type + 1)
        self._type_names = [''] * (n_atom_type + 1)

        for i, (iline, section) in enumerate(_iline_section):
            if i != len(_iline_section) - 1:
                iline_next = _iline_section[i + 1][0]
            else:
                iline_next = len(lines)
            if section == 'Masses':
                self.parse_masses(lines[iline + 2: iline_next], n_atom_type)  # there is a blank line after each section header
            if section == 'Atoms':
                self.parse_atoms(lines[iline + 2: iline_next], n_atom)
            if section == 'Bonds':
                self.parse_bonds(lines[iline + 2: iline_next], n_bond)

    def parse_header(self, lines):
        n_atom = n_bond = n_angle = n_dihedral = n_improper = n_atom_type = 0
        for line in lines:
            if line.endswith(' atoms'):
                n_atom = int(line.split()[0])
            elif line.endswith(' bonds'):
                n_bond = int(line.split()[0])
            elif line.endswith(' angles'):
                n_angle = int(line.split()[0])
            elif line.endswith(' dihedrals'):
                n_dihedral = int(line.split()[0])
            elif line.endswith(' impropers'):
                n_improper = int(line.split()[0])
            elif line.endswith(' atom types'):
                n_atom_type = int(line.split()[0])
            # elif line.endswith(' bond types'):
            #     n_bond_type = int(line.split()[0])
            # elif line.endswith(' angle types'):
            #     n_angle_type = int(line.split()[0])
            # elif line.endswith(' dihedral types'):
            #     n_dihedral_type = int(line.split()[0])
            # elif line.endswith(' improper types'):
            #     n_improper_type = int(line.split()[0])
            else:
                continue

        return n_atom, n_bond, n_angle, n_dihedral, n_improper, n_atom_type

    def parse_masses(self, lines, n_atom_type):
        for i in range(n_atom_type):
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

    def parse_atoms(self, lines, n_atom):
        atoms = [Atom() for i in range(n_atom)]
        # i don't know how many molecules there will be
        # but it can not be larger than n_atom
        # i'll kick out the redundant molecules later
        molecules = [Molecule() for i in range(n_atom)]
        for i in range(n_atom):
            words = lines[i].strip().split()
            atom_id = int(words[0])
            mol_id = int(words[1])
            type_id = int(words[2])
            charge = float(words[3])

            mol = molecules[mol_id-1]
            if not mol.initiated:
                mol.id = mol_id - 1  # mol.id starts from 0
                mol.name = words[9] if (len(words) > 7 and words[7] == '#') else 'UNK'

            atom = atoms[atom_id - 1]
            atom.id = atom_id - 1  # atom.id starts from 0
            mol.add_atom(atom)
            atom.charge = charge
            atom.mass = self._type_masses[type_id]
            atom.type = self._type_names[type_id]
            atom.symbol = Element.guess_from_atom_type(atom.type).symbol


        # kick out redundant molecules
        molecules = [mol for mol in molecules if mol.initiated]
        for mol in molecules:
            for i, atom in enumerate(mol.atoms):
                atom.name = atom.symbol + str(i + 1)  # atomic symbol + index inside mol starting from 1

        self.init_from_molecules(molecules)

    def parse_bonds(self, lines, n_bond):
        pass

    def parse_angles(self, lines, n_angle):
        pass

    def parse_dihedrals(self, lines, n_dihedral):
        pass

    def parse_impropers(self, lines, n_improper):
        pass


