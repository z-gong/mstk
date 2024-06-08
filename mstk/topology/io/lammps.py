import numpy as np
from mstk.topology.atom import Atom
from mstk.topology.molecule import Molecule
from mstk.topology.topology import Topology
from mstk.chem.element import Element
from mstk import logger


class LammpsData():
    '''
    Generate topology from Lammps data file.

    Atom masses, types, charges, connectivity, unit cell and atom positions are parsed.
    The velocities and force field parameters are ignored.

    The string of atom types are assumed to be the comment after the mass.
    If such comment is not provided, then the string of atom types will be 'AT1', 'AT2'...
    If any atom type has mass smaller than 1.0 and atom type comment started with 'D',
    then it will be treated as Drude particle. e.g.

    >>> Masses
    >>>   1 1.008 # HC
    >>>   2 16.00 # OH
    >>>   3 0.800 # D

    HC and OH will be treated as the atom type string of the first and second atom types.
    The third atom type will be treated as Drude particle.

    Care should be taken for the impropers.
    The central atom for impropers in mstk is the first atom.
    But in Lammps data file, the central atom can be the first, second or third atom,
    depending on the `improper_style`.
    The argument `improper_center` controls which atom should be treated as the center atom,
    which is set to 1 by default.

    Parameters
    ----------
    file : str
    improper_center : [1, 2, 3, 4]
    kwargs : dict
    Ignored

    Examples
    --------

    >>> lammps = LammpsData('data.lmp')
    >>> topology = lammps.topology

    For `improper_style cvff`, the central atom is the third.
    Then the LammpsData should be loaded as

    >>> lammps = LammpsData('data.lmp', improper_center=3)

    Attributes
    ----------
    topology : Topology
        Parsed topology from the data
    n_atom_type : int
        Number of atom types
    type_names : list of str
        Names of these atom types
    type_masses : list of float
        Masses of these atom types
    xlo : float
        The origin of the cell in x direction. Not really useful
    ylo : float
        The origin of the cell in z direction. Not really useful
    zlo : float
        The origin of the cell in z direction. Not really useful
    '''

    def __init__(self, file, improper_center=None, **kwargs):
        self.topology = Topology()
        self._parse(file, improper_center)

    def _parse(self, file, improper_center):
        with open(file) as f:
            lines = f.read().splitlines()

        self.topology.remark = lines[0].strip()

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

        n_atom, n_bond, n_angle, n_dihedral, n_improper = \
            self._parse_header(lines[:_iline_section[0][0]])

        self.type_masses = [0.] * (self.n_atom_type + 1)
        self.type_names = [''] * (self.n_atom_type + 1)

        for i, (iline, section) in enumerate(_iline_section):
            if i != len(_iline_section) - 1:
                iline_next = _iline_section[i + 1][0]
            else:
                iline_next = len(lines)
            if section == 'Masses':
                # there is a blank line after each section header
                self._parse_atom_types(lines[iline + 2: iline_next])
            if section == 'Atoms':
                self._parse_atoms(lines[iline + 2: iline_next], n_atom)
            if section == 'Bonds':
                self._parse_bonds(lines[iline + 2: iline_next], n_bond)
            if section == 'Angles':
                self._parse_angles(lines[iline + 2: iline_next], n_angle)
            if section == 'Dihedrals':
                self._parse_dihedrals(lines[iline + 2: iline_next], n_dihedral)
            if section == 'Impropers':
                self._parse_impropers(lines[iline + 2: iline_next], n_improper, improper_center)

    def _parse_header(self, lines):
        n_atom = n_bond = n_angle = n_dihedral = n_improper = self.n_atom_type = 0
        la = lb = lc = bx = cx = cy = 0
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
                self.n_atom_type = int(words[0])
            elif line.endswith(' xlo xhi'):
                la = (float(words[1]) - float(words[0])) / 10
                self.xlo = float(words[0]) / 10
            elif line.endswith(' ylo yhi'):
                lb = (float(words[1]) - float(words[0])) / 10
                self.ylo = float(words[0]) / 10
            elif line.endswith(' zlo zhi'):
                lc = (float(words[1]) - float(words[0])) / 10
                self.zlo = float(words[0]) / 10
            elif line.endswith(' xy xz yz'):
                bx, cx, cy = tuple(map(lambda x: float(x) / 10, words[:3]))
            else:
                continue

        self.topology.cell.set_box([[la, 0, 0], [bx, lb, 0], [cx, cy, lc]])

        return n_atom, n_bond, n_angle, n_dihedral, n_improper

    def _parse_atom_types(self, lines):
        for i in range(self.n_atom_type):
            words = lines[i].strip().split()
            type_id = int(words[0])
            mass = float(words[1])
            self.type_masses[type_id] = mass
            # TODO This is not robust but there's no better way
            if len(words) > 3 and words[2] == '#':
                if words[-1].startswith('D') and mass < 1:
                    self.type_names[type_id] = 'DP_'
                    logger.warning(f'Atom type {type_id} is considered to be Drude particle')
                else:
                    self.type_names[type_id] = words[3]
            else:
                self.type_names[type_id] = 'AT' + str(type_id)

    def _parse_atoms(self, lines, n_atom):
        atoms = [Atom() for i in range(n_atom)]
        mol_names: {int: str} = {}
        box = self.topology.cell.get_size()
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

            if mol_id not in mol_names:
                mol_names[mol_id] = words[9] if (
                        len(words) > 9 and words[7] == '#') else 'M%i' % mol_id

            atom = atoms[atom_id - 1]
            atom._mol_id = mol_id  # temporary attribute for identifying molecules
            atom.id = atom_id - 1  # atom.id starts from 0
            atom.charge = charge
            atom.mass = self.type_masses[type_id]
            atom.type = self.type_names[type_id]
            if atom.type == 'DP_':
                atom.is_drude = True
                atom.symbol = 'DP'
            else:
                atom.symbol = Element.guess_from_atom_type(atom.type).symbol
            atom.position = np.array([x - self.xlo + ix * box[0],
                                      y - self.ylo + iy * box[1],
                                      z - self.zlo + iz * box[2]])

        molecules = [Molecule(name) for id, name in sorted(mol_names.items())]
        for atom in sorted(atoms, key=lambda x: x.id):
            mol = molecules[atom._mol_id - 1]
            mol.add_atom(atom, update_topology=False)
        for mol in molecules:
            for i, atom in enumerate(mol.atoms):
                # atomic symbol + index inside mol starting from 1
                atom.name = atom.symbol + str(i + 1)

        self.topology.update_molecules(molecules)

    def _parse_bonds(self, lines, n_bond):
        atoms = self.topology.atoms
        for i in range(n_bond):
            words = lines[i].strip().split()
            atom1 = atoms[int(words[2]) - 1]
            atom2 = atoms[int(words[3]) - 1]
            atom1.molecule.add_bond(atom1, atom2)

    def _parse_angles(self, lines, n_angle):
        atoms = self.topology.atoms
        for i in range(n_angle):
            words = lines[i].strip().split()
            atom1 = atoms[int(words[2]) - 1]
            atom2 = atoms[int(words[3]) - 1]
            atom3 = atoms[int(words[4]) - 1]
            atom1.molecule.add_angle(atom1, atom2, atom3)

    def _parse_dihedrals(self, lines, n_dihedral):
        atoms = self.topology.atoms
        for i in range(n_dihedral):
            words = lines[i].strip().split()
            atom1 = atoms[int(words[2]) - 1]
            atom2 = atoms[int(words[3]) - 1]
            atom3 = atoms[int(words[4]) - 1]
            atom4 = atoms[int(words[5]) - 1]
            atom1.molecule.add_dihedral(atom1, atom2, atom3, atom4)

    def _parse_impropers(self, lines, n_improper, improper_center):
        if improper_center is None:
            improper_center = 1
            logger.warning('improper_center undefined, '
                           'will treat the first atom as the central atom for impropers')
        elif improper_center not in (1, 2, 3, 4):
            raise Exception('improper_center should be a integer from 1 to 4')

        atoms = self.topology.atoms
        for i in range(n_improper):
            words = lines[i].strip().split()
            sides = []
            sides.append(atoms[int(words[2]) - 1])
            sides.append(atoms[int(words[3]) - 1])
            sides.append(atoms[int(words[4]) - 1])
            sides.append(atoms[int(words[5]) - 1])
            center = sides.pop(improper_center - 1)
            sides[0].molecule.add_improper(center, *sides)


Topology.register_format('.lmp', LammpsData)
