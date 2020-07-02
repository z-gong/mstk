import numpy as np
from .atom import Atom
from .molecule import Molecule
from .topology import Topology
from ..forcefield import Element
from .. import logger


class LammpsData():
    '''
    Generate Topology from Lammps data file.

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
    The central atom for impropers in mstools is the first atom.
    But in Lammps data file, the central atom can be the first, second or third atom,
    depending on the `improper_style`.
    The argument `improper_center` controls which atom should be treated as the center atom,
    which is set to 1 by default.
    e.g. for `improper_style cvff`, the central atom is the third.
    Then the LammpsData should be loaded as

    >>> topology = LammpsData.read('data.lmp', improper_center=3)

    '''

    @staticmethod
    def read(file, improper_center=None, **kwargs):
        '''
        Parse a Lammps data file

        Parameters
        ----------
        file : str
        improper_center : [1, 2, 3, 4]
        kwargs : dict
            Ignored
        '''
        with open(file) as f:
            lines = f.read().splitlines()

        topology = Topology()
        topology.remark = lines[0].strip()

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

        n_atom, n_bond, n_angle, n_dihedral, n_improper, n_atom_type, xlo, ylo, zlo = \
            LammpsData._parse_header(topology, lines[:_iline_section[0][0]])

        _type_masses = [0.] * (n_atom_type + 1)
        _type_names = [''] * (n_atom_type + 1)

        for i, (iline, section) in enumerate(_iline_section):
            if i != len(_iline_section) - 1:
                iline_next = _iline_section[i + 1][0]
            else:
                iline_next = len(lines)
            if section == 'Masses':
                # there is a blank line after each section header
                LammpsData._parse_atom_types(lines[iline + 2: iline_next], n_atom_type,
                                             _type_names, _type_masses)
            if section == 'Atoms':
                LammpsData._parse_atoms(topology, lines[iline + 2: iline_next], n_atom,
                                        _type_names, _type_masses, xlo, ylo, zlo)
            if section == 'Bonds':
                LammpsData._parse_bonds(topology, lines[iline + 2: iline_next], n_bond)
            if section == 'Angles':
                LammpsData._parse_angles(topology, lines[iline + 2: iline_next], n_angle)
            if section == 'Dihedrals':
                LammpsData._parse_dihedrals(topology, lines[iline + 2: iline_next], n_dihedral)
            if section == 'Impropers':
                LammpsData._parse_impropers(topology, lines[iline + 2: iline_next], n_improper,
                                            improper_center)

        return topology

    @staticmethod
    def _parse_header(top, lines):
        n_atom = n_bond = n_angle = n_dihedral = n_improper = n_atom_type = 0
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
                n_atom_type = int(words[0])
            elif line.endswith(' xlo xhi'):
                la = (float(words[1]) - float(words[0])) / 10
                xlo = float(words[0]) / 10
            elif line.endswith(' ylo yhi'):
                lb = (float(words[1]) - float(words[0])) / 10
                ylo = float(words[0]) / 10
            elif line.endswith(' zlo zhi'):
                lc = (float(words[1]) - float(words[0])) / 10
                zlo = float(words[0]) / 10
            elif line.endswith(' xy xz yz'):
                bx, cx, cy = tuple(map(lambda x: float(x) / 10, words[:3]))
            else:
                continue

        top.cell.set_box([[la, 0, 0], [bx, lb, 0], [cx, cy, lc]])

        return n_atom, n_bond, n_angle, n_dihedral, n_improper, n_atom_type, xlo, ylo, zlo

    @staticmethod
    def _parse_atom_types(lines, n_atom_type, _type_names, _type_masses):
        for i in range(n_atom_type):
            words = lines[i].strip().split()
            type_id = int(words[0])
            mass = float(words[1])
            _type_masses[type_id] = mass
            # TODO This is not robust but there's no better way
            if len(words) > 3 and words[2] == '#':
                if words[-1].startswith('D') and mass < 1:
                    _type_names[type_id] = 'DP_'
                    logger.warning(f'Atom type {type_id} is considered to be Drude particle')
                else:
                    _type_names[type_id] = words[3]
            else:
                _type_names[type_id] = 'AT' + str(type_id)

    @staticmethod
    def _parse_atoms(top, lines, n_atom, _type_names, _type_masses, _xlo, _ylo, _zlo):
        atoms = [Atom() for i in range(n_atom)]
        mol_names: {int: str} = {}
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
            atom.mass = _type_masses[type_id]
            atom.type = _type_names[type_id]
            if atom.type == 'DP_':
                atom.is_drude = True
                atom.symbol = 'DP'
            else:
                atom.symbol = Element.guess_from_atom_type(atom.type).symbol
            atom.position = np.array([x - _xlo + ix * top.cell.size[0],
                                      y - _ylo + iy * top.cell.size[1],
                                      z - _zlo + iz * top.cell.size[2]])

        molecules = [Molecule(name) for id, name in sorted(mol_names.items())]
        for atom in sorted(atoms, key=lambda x: x.id):
            mol = molecules[atom._mol_id - 1]
            mol.add_atom(atom)
        for mol in molecules:
            for i, atom in enumerate(mol.atoms):
                # atomic symbol + index inside mol starting from 1
                atom.name = atom.symbol + str(i + 1)

        top.update_molecules(molecules)

    @staticmethod
    def _parse_bonds(top, lines, n_bond):
        atoms = top.atoms
        for i in range(n_bond):
            words = lines[i].strip().split()
            atom1 = atoms[int(words[2]) - 1]
            atom2 = atoms[int(words[3]) - 1]
            atom1.molecule.add_bond(atom1, atom2)

    @staticmethod
    def _parse_angles(top, lines, n_angle):
        atoms = top.atoms
        for i in range(n_angle):
            words = lines[i].strip().split()
            atom1 = atoms[int(words[2]) - 1]
            atom2 = atoms[int(words[3]) - 1]
            atom3 = atoms[int(words[4]) - 1]
            atom1.molecule.add_angle(atom1, atom2, atom3)

    @staticmethod
    def _parse_dihedrals(top, lines, n_dihedral):
        atoms = top.atoms
        for i in range(n_dihedral):
            words = lines[i].strip().split()
            atom1 = atoms[int(words[2]) - 1]
            atom2 = atoms[int(words[3]) - 1]
            atom3 = atoms[int(words[4]) - 1]
            atom4 = atoms[int(words[5]) - 1]
            atom1.molecule.add_dihedral(atom1, atom2, atom3, atom4)

    @staticmethod
    def _parse_impropers(top, lines, n_improper, improper_center):
        if improper_center is None:
            improper_center = 1
            logger.warning('improper_center undefined, '
                           'will treat the first atom as the central atom for impropers')
        elif improper_center not in (1, 2, 3, 4):
            raise Exception('improper_center should be a integer from 1 to 4')

        atoms = top.atoms
        for i in range(n_improper):
            words = lines[i].strip().split()
            sides = []
            sides.append(atoms[int(words[2]) - 1])
            sides.append(atoms[int(words[3]) - 1])
            sides.append(atoms[int(words[4]) - 1])
            sides.append(atoms[int(words[5]) - 1])
            center = sides.pop(improper_center - 1)
            sides[0].molecule.add_improper(center, *sides)
