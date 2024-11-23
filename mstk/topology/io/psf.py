import math
from mstk.chem.element import Element
from mstk.topology.atom import Atom
from mstk.topology.molecule import Molecule
from mstk.topology.topology import Topology
from mstk.topology.virtualsite import TIP4PSite


class Psf:
    '''
    Generate Topology from PSF file.

    Atoms, bonds, angles, dihedrals, impropers and virtual sites are parsed.
    Polarizability and thole parameters for Drude model are parsed,
    but anisotroic polarizability are ignored.
    All others like hydrogen bonds are ignored.

    Parameters
    ----------
    file : str
    split_molecule : str
        Can be 'residue', 'whole', 'atom'. Default is 'residue'.
        If set to 'residue', the topology will be split into molecules based on the bonds between residues.
        If set to 'whole', all the atoms will be put into one molecule.
        If set to 'atom', the topology will be split into molecules based on the bonds between atoms.

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
        self._molecule = Molecule()  # a single molecule holds all the atoms
        self._parse(file, **kwargs)

    def _parse(self, file, split_molecule='residue'):
        with open(file) as f:
            lines = f.read().splitlines()

        words = lines[0].strip().split()
        if 'PSF' not in words:
            raise Exception('Invalid psf file, PSF not found in first line')
        is_drude = 'DRUDE' in words

        # there are other sections, but i don't care
        _sections = ['NTITLE', 'NATOM', 'NBOND', 'NTHETA', 'NPHI', 'NIMPHI', 'NUMLP NUMLPH']

        # record the line number (starts from 0) of section headers
        _iline_section = []
        for i, line in enumerate(lines):
            words = line.strip().split('!')
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
            if section == 'NTITLE':
                # the number at section line indicates the line counts of the section
                self._parse_title(lines[iline: iline_next])
            if section == 'NATOM':
                self._parse_atoms(lines[iline: iline_next], is_drude)
            if section == 'NBOND':
                self._parse_bonds(lines[iline: iline_next])
            if section == 'NTHETA':
                self._parse_angles(lines[iline: iline_next])
            if section == 'NPHI':
                self._parse_dihedrals(lines[iline: iline_next])
            if section == 'NIMPHI':
                self._parse_impropers(lines[iline: iline_next])
            if section == 'NUMLP NUMLPH':
                self._parse_virtual_sites(lines[iline: iline_next])

        if split_molecule == 'residue':
            self.topology.update_molecules(self._molecule.split_residues(consecutive=True))
        elif split_molecule == 'whole':
            self.topology.update_molecules([self._molecule])
        elif split_molecule == 'atom':
            self.topology.update_molecules(self._molecule.split(consecutive=True))
        else:
            raise Exception(f'Invalid option for split_molecule: {split_molecule}')

    def _parse_title(self, lines):
        nline = int(lines[0].strip().split()[0])
        self.topology.remark = '\n'.join(lines[1:1 + nline]).strip()

    def _parse_atoms(self, lines, parse_drude):
        '''
        in psf file, drude particles always appear after attached atoms
        '''
        mol = self._molecule
        n_atom = int(lines[0].strip().split()[0])
        res_names = {}  # {int: str}
        res_atoms = {}  # {int: [Atom]}
        for i in range(n_atom):
            words = lines[1 + i].strip().split()
            res_id = int(words[2])
            res_name = words[3]
            atom_name = words[4]
            atom_type = words[5]
            charge, mass = list(map(float, words[6:8]))
            if parse_drude:
                alpha, thole = list(map(float, words[9:11]))
            else:
                alpha = thole = 0.0

            atom = Atom(atom_name)
            atom.charge = charge
            atom.mass = mass
            atom.type = atom_type
            atom.alpha = -alpha / 1000  # convert A^3 to nm^3
            atom.thole = thole * 2
            if parse_drude and atom.type.startswith('D'):
                atom.is_drude = True
                atom.symbol = 'DP'
            else:
                atom.symbol = Element.guess_from_atom_type(atom.type).symbol

            mol.add_atom(atom)

            if res_id not in res_names:
                res_names[res_id] = res_name
                res_atoms[res_id] = []
            res_atoms[res_id].append(atom)

        for resid, resname in res_names.items():
            atoms = res_atoms[resid]
            mol.add_residue(resname, atoms, refresh_residues=False)
        mol.refresh_residues()

    def _parse_bonds(self, lines):
        n_bond = int(lines[0].strip().split()[0])
        atoms = self._molecule.atoms
        for i in range(math.ceil(n_bond / 4)):
            words = lines[i + 1].strip().split()
            for j in range(4):
                if i * 4 + j + 1 > n_bond:
                    break
                atom1, atom2 = atoms[int(words[j * 2]) - 1], atoms[int(words[j * 2 + 1]) - 1]
                atom1.molecule.add_bond(atom1, atom2)

    def _parse_angles(self, lines):
        n_angle = int(lines[0].strip().split()[0])
        atoms = self._molecule.atoms
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
        atoms = self._molecule.atoms
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
        atoms = self._molecule.atoms
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

    def _parse_virtual_sites(self, lines):
        '''
        Currently only TIP4P style virtual sites are parsed
        '''
        n_vsite, n_holder = map(int, lines[0].strip().split()[:2])
        if n_holder != 4 * n_vsite:
            raise Exception('Only TIP4P style virtual sites supported: ' + lines[0])

        d_list = [0] * n_vsite
        for i in range(n_vsite):
            words = lines[i + 1].strip().split()
            if len(words) != 6 or words[0] != '3' or words[2] != 'F':
                raise Exception('Only TIP4P style virtual sites supported: ' + lines[i + 1])
            d, theta, phi = map(float, words[-3:])
            if theta != 0 or phi != 0:
                raise Exception('Only TIP4P style virtual sites supported: ' + lines[i + 1])
            d_list[i] = -d / 10  # convert from A to nm

        atoms = self._molecule.atoms
        for i in range(math.ceil(n_vsite / 2)):
            words = lines[i + 1 + n_vsite].strip().split()
            for j in range(2):
                if i * 2 + j + 1 > n_vsite:
                    break
                atom1, atom2, atom3, atom4 = atoms[int(words[j * 4 + 0]) - 1], \
                                             atoms[int(words[j * 4 + 1]) - 1], \
                                             atoms[int(words[j * 4 + 2]) - 1], \
                                             atoms[int(words[j * 4 + 3]) - 1]
                atom1.virtual_site = TIP4PSite([atom2, atom3, atom4], [d_list[i * 2 + j]])

    @staticmethod
    def save_to(top, file, **kwargs):
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

        string = ''

        if top.is_drude:
            string += 'PSF DRUDE\n\n'
        else:
            string += 'PSF\n\n'

        string += '%8i !NTITLE\n' % 1
        string += ' REMARKS Created by mstk\n\n'

        string += '%8i !NATOM\n' % top.n_atom
        for i, atom in enumerate(top.atoms):
            string += '%8i %-4s %-4i %-7s %-4s %-4s %10.6f %13.4f %11i %8.4f %8.4f\n' % (
                atom.id + 1, 'S', atom.residue.id + 1, atom.residue.name, atom.name, atom.type,
                atom.charge, atom.mass, 0, -atom.alpha * 1000, atom.thole / 2
            )
        string += '\n'

        n_bond = top.n_bond
        string += '%8i !NBOND: bonds\n' % n_bond
        for i, bond in enumerate(top.bonds):
            string += '%8i%8i' % (bond.atom1.id + 1, bond.atom2.id + 1)
            if (i + 1) % 4 == 0 or i + 1 == n_bond:
                string += '\n'
        string += '\n'

        n_angle = top.n_angle
        string += '%8i !NTHETA: angles\n' % n_angle
        for i, angle in enumerate(top.angles):
            string += '%8i%8i%8i' % (angle.atom1.id + 1, angle.atom2.id + 1,
                                     angle.atom3.id + 1)
            if (i + 1) % 3 == 0 or i + 1 == n_angle:
                string += '\n'
        string += '\n'

        n_dihedral = top.n_dihedral
        string += '%8i !NPHI: dihedrals\n' % n_dihedral
        for i, dihedral in enumerate(top.dihedrals):
            string += '%8i%8i%8i%8i' % (dihedral.atom1.id + 1, dihedral.atom2.id + 1,
                                        dihedral.atom3.id + 1, dihedral.atom4.id + 1)
            if (i + 1) % 2 == 0 or i + 1 == n_dihedral:
                string += '\n'
        string += '\n'

        n_improper = top.n_improper
        string += '%8i !NIMPHI: impropers\n' % n_improper
        for i, improper in enumerate(top.impropers):
            string += '%8i%8i%8i%8i' % (improper.atom1.id + 1, improper.atom2.id + 1,
                                        improper.atom3.id + 1, improper.atom4.id + 1)
            if (i + 1) % 2 == 0 or i + 1 == n_improper:
                string += '\n'
        string += '\n'

        string += '%8i !NDON: donors\n\n' % 0
        string += '%8i !NACC: acceptors\n\n' % 0
        string += '%8i !NNB\n\n\n' % 0

        vsite_pairs = top.get_virtual_site_pairs()
        n_vsite = len(vsite_pairs)
        string += '%8i %i !NUMLP NUMLPH\n' % (n_vsite, 4 * n_vsite)
        for parent, avsite in vsite_pairs:
            if type(avsite.virtual_site) is not TIP4PSite:
                raise Exception('Virtual sites other than TIP4PSite haven\'t been implemented')
            string += '%8i%8i  F %10.6f %10.6f %10.6f\n' % (
                3, parent.id + 1, -avsite.virtual_site.parameters[0] * 10, 0, 0)
        for i, (parent, avsite) in enumerate(vsite_pairs):
            O, H1, H2 = avsite.virtual_site.parents
            string += '%8i%8i%8i%8i' % (avsite.id + 1, O.id + 1, H1.id + 1, H2.id + 1)
            if (i + 1) % 2 == 0 or i + 1 == n_vsite:
                string += '\n'
        string += '\n'

        string += '%8i !NUMANISO\n\n' % 0

        with open(file, 'wb') as f:
            f.write(string.encode())


Topology.register_format('.psf', Psf)
