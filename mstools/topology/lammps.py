import numpy as np
from ..forcefield import Element
from ..topology import Atom, Molecule, Topology


class LammpsData(Topology):
    '''
    Parse the topology information from Lammps data file
    The box and atom positions are also parsed,
    the velocities and force field parameters are ignored
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

        n_atom, n_bond, n_angle, n_dihedral, n_improper, n_atom_type = self.parse_header(
            lines[:_iline_section[0][0]])

        self._type_masses = [0.] * (n_atom_type + 1)
        self._type_names = [''] * (n_atom_type + 1)

        for i, (iline, section) in enumerate(_iline_section):
            if i != len(_iline_section) - 1:
                iline_next = _iline_section[i + 1][0]
            else:
                iline_next = len(lines)
            if section == 'Masses':
                # there is a blank line after each section header
                self.parse_masses(lines[iline + 2: iline_next], n_atom_type)
            if section == 'Atoms':
                self.parse_atoms(lines[iline + 2: iline_next], n_atom)
            if section == 'Bonds':
                self.parse_bonds(lines[iline + 2: iline_next], n_bond)
            if section == 'Angles':
                self.parse_angles(lines[iline + 2: iline_next], n_angle)
            if section == 'Dihedrals':
                self.parse_dihedrals(lines[iline + 2: iline_next], n_dihedral)
            if section == 'Impropers':
                self.parse_impropers(lines[iline + 2: iline_next], n_improper)

    def parse_header(self, lines):
        n_atom = n_bond = n_angle = n_dihedral = n_improper = n_atom_type = 0
        box = [0., 0., 0.]
        for line in lines:
            line = line.strip()
            words = line.split()
            if line.endswith(' atoms'):
                n_atom = int(words[0])
            elif line.endswith(' bonds'):
                n_bond = int(words[0])
            elif line.endswith(' angles'):
                n_angle = int(words[0])
            elif line.endswith(' dihedrals'):
                n_dihedral = int(words[0])
            elif line.endswith(' impropers'):
                n_improper = int(words[0])
            elif line.endswith(' atom types'):
                n_atom_type = int(words[0])
            elif line.endswith(' xlo xhi'):
                box[0] = (float(words[1]) - float(words[0])) / 10
                self._xlo = float(words[0]) / 10
            elif line.endswith(' ylo yhi'):
                box[1] = (float(words[1]) - float(words[0])) / 10
                self._ylo = float(words[0]) / 10
            elif line.endswith(' zlo zhi'):
                box[2] = (float(words[1]) - float(words[0])) / 10
                self._zlo = float(words[0]) / 10
            else:
                continue

        self.box = box

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
            x = float(words[4]) / 10  # convert from A to nm
            y = float(words[5]) / 10
            z = float(words[6]) / 10
            try:
                ix = int(words[7])
                iy = int(words[8])
                iz = int(words[9])
            except:
                ix = iy = iz = 0

            mol = molecules[mol_id - 1]
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
            atom.has_position = True
            atom.position = np.array([x - self._xlo + ix * self.box[0],
                                      y - self._ylo + iy * self.box[1],
                                      z - self._zlo + iz * self.box[2]])

        # kick out redundant molecules
        molecules = [mol for mol in molecules if mol.initiated]
        for mol in molecules:
            for i, atom in enumerate(mol.atoms):
                # atomic symbol + index inside mol starting from 1
                atom.name = atom.symbol + str(i + 1)

        self.init_from_molecules(molecules)

    def parse_bonds(self, lines, n_bond):
        atoms = self.atoms
        for i in range(n_bond):
            words = lines[i].strip().split()
            atom1 = atoms[int(words[2]) - 1]
            atom2 = atoms[int(words[3]) - 1]
            atom1.molecule.add_bond(atom1, atom2)

    def parse_angles(self, lines, n_angle):
        atoms = self.atoms
        for i in range(n_angle):
            words = lines[i].strip().split()
            atom1 = atoms[int(words[2]) - 1]
            atom2 = atoms[int(words[3]) - 1]
            atom3 = atoms[int(words[4]) - 1]
            atom1.molecule.add_angle(atom1, atom2, atom3)

    def parse_dihedrals(self, lines, n_dihedral):
        atoms = self.atoms
        for i in range(n_dihedral):
            words = lines[i].strip().split()
            atom1 = atoms[int(words[2]) - 1]
            atom2 = atoms[int(words[3]) - 1]
            atom3 = atoms[int(words[4]) - 1]
            atom4 = atoms[int(words[5]) - 1]
            atom1.molecule.add_dihedral(atom1, atom2, atom3, atom4)

    def parse_impropers(self, lines, n_improper):
        atoms = self.atoms
        for i in range(n_improper):
            words = lines[i].strip().split()
            atom1 = atoms[int(words[2]) - 1]
            atom2 = atoms[int(words[3]) - 1]
            atom3 = atoms[int(words[4]) - 1]
            atom4 = atoms[int(words[5]) - 1]
            atom1.molecule.add_improper(atom1, atom2, atom3, atom4)
