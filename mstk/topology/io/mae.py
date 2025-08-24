import gzip
import shlex
from mstk.topology import Atom, Bond, Molecule, Topology, UnitCell
from mstk.chem.element import Element
from mstk.chem.constant import DEG2RAD, RAD2DEG
from mstk.utils import flatten


class MaeBlock:
    def __init__(self, name):
        self.name = name
        self.data = {}
        self.parent = None
        self.children = []

    def add_child(self, child):
        self.children.append(child)
        child.parent = self

    def __repr__(self):
        return f'<MaeBlock: {self.name}>'


class Mae:
    '''
    Generate Topology from an atomistic MAE or MAEGZ file of Schrödinger.

    The atom symbols, positions, bonds are parsed.
    Atom types, charges, formal charges, residues will be parsed if avaiable.
    Atom names will be parsed only if atom names for all atoms in one entry are available. Otherwise, atoms names in this entry will be re-assigned.
    The first available unit cell from the structure entryies will be parsed.
    All other information are ignored.

    An MAE file can contain multiple structures. Each structure entry can contain multiple molecules.
    The parser will split each structure entry into molecules based on the bonds between residues or atoms.
    The topology will contain all the molecules in all the structure entries.

    Parameters
    ----------
    file : str
        MAE file
    split_molecule : str
        Can be 'residue', 'whole', 'atom'. Default is 'residue'.
        If set to 'residue', each structure entry will be split into molecules based on the bonds between residues.
        If set to 'whole', each structure entry will be put into one molecule.
        If set to 'atom', each structure entry will be split into molecules based on the bonds between atoms.


    Attributes
    ----------
    topology : Topology

    Examples
    --------
    >>> mae = Mae('input.mae')
    >>> topology = mae.topology

    '''

    def __init__(self, file, **kwargs):
        suffix = file.split('.')[-1].lower()
        if suffix == 'mae':
            f = open(file)
        elif suffix == 'maegz':
            f = gzip.open(file, 'rt')
        else:
            raise Exception('Input must be a MAE or MAEGZ file')
        self.lines = f.readlines()
        f.close()

        self.root = MaeBlock('root')
        self.current_line = 0
        while self.current_line < len(self.lines):
            self._parse(self.root)

        self.topology = self._build_topology(**kwargs)

    def _parse(self, block):
        while self.current_line < len(self.lines):
            line = self.lines[self.current_line].strip()
            self.current_line += 1
            if line == '{':
                self._parse_property_block(block)
            if line == 'f_m_ct {':
                block = MaeBlock('f_m_ct')
                self.root.add_child(block)
                self._parse_property_block(block)

    def _parse_property_block(self, block):
        attributes = []
        values = []
        _value = False
        while self.current_line < len(self.lines):
            line = self.lines[self.current_line].strip()
            self.current_line += 1
            if line == '}':
                block.data = {a: v for a, v in zip(attributes, values)}
                break
            if line.endswith('{'):
                name = line.split()[0]
                if '[' in name:
                    name = name.split('[')[0]
                    parse = self._parse_list_block
                else:
                    parse = self._parse_property_block
                child = MaeBlock(name)
                block.add_child(child)
                block.data = {a: v for a, v in zip(attributes, values)}
                parse(child)
            if line == ':::':
                _value = True
                continue
            if _value:
                values.append(line.strip('"'))
            else:
                attributes.append(line)

    def _parse_list_block(self, block):
        attributes = []
        rows = []
        _row = False
        while self.current_line < len(self.lines):
            line = self.lines[self.current_line].strip()
            self.current_line += 1
            if line == '}':
                data = {a: [] for a in attributes}
                for row in rows:
                    for j, attr in enumerate(attributes):
                        data[attr].append(row[j])
                block.data = data
                break
            if line == ':::':
                _row = True
                continue
            if _row:
                rows.append(shlex.split(line))
            else:
                attributes.append(line)

    def _build_topology(self, split_molecule='residue'):
        molecules = []
        cell = UnitCell()
        for f_m_ct in self.root.children:
            if cell.volume == 0 and 'r_pdb_PDB_CRYST1_a' in f_m_ct.data:
                a, b, c, alpha, beta, gamma = float(f_m_ct.data['r_pdb_PDB_CRYST1_a']) / 10, \
                                              float(f_m_ct.data['r_pdb_PDB_CRYST1_b']) / 10, \
                                              float(f_m_ct.data['r_pdb_PDB_CRYST1_c']) / 10, \
                                              float(f_m_ct.data['r_pdb_PDB_CRYST1_alpha']) * DEG2RAD, \
                                              float(f_m_ct.data['r_pdb_PDB_CRYST1_beta']) * DEG2RAD, \
                                              float(f_m_ct.data['r_pdb_PDB_CRYST1_gamma']) * DEG2RAD
                cell = UnitCell([[a, b, c], [alpha, beta, gamma]])

            atom_block = next(block for block in f_m_ct.children if block.name == 'm_atom')
            residue_numbers = atom_block.data.get('i_m_residue_number', None)  # optional
            residue_names = atom_block.data.get('s_m_pdb_residue_name', None)  # optional
            formal_charges = atom_block.data.get('i_m_formal_charge', None)  # optional
            atom_names = atom_block.data.get('s_m_atom_name', None)  # optional
            atom_types = atom_block.data.get('s_m_atom_type', None)  # optional
            charges = atom_block.data.get('r_m_charge', None)  # optional

            mol = Molecule(name=f_m_ct.data.get('s_m_title', 'UNK'))
            for i, atom_number in enumerate(atom_block.data['i_m_atomic_number']):
                element = Element(int(atom_number))
                atom = Atom()
                atom.symbol = element.symbol
                atom.mass = element.mass
                mol.add_atom(atom)

                x, y, z = float(atom_block.data['r_m_x_coord'][i]) / 10, \
                          float(atom_block.data['r_m_y_coord'][i]) / 10, \
                          float(atom_block.data['r_m_z_coord'][i]) / 10
                atom.position = [x, y, z]

                if formal_charges:
                    atom.formal_charge = int(formal_charges[i])
                if atom_types:
                    atom.type = atom_types[i]
                if charges:
                    atom.charge = float(charges[i])

            if residue_numbers:
                groups = []
                current_group = [0]
                for i in range(1, mol.n_atom):
                    if residue_numbers[i] == residue_numbers[i - 1]:
                        current_group.append(i)
                    else:
                        groups.append(current_group)
                        current_group = [i]
                groups.append(current_group)

                for group in groups:
                    atoms = [mol.atoms[i] for i in group]
                    resname = residue_names[group[0]].strip()  # residue name in MAE are filled with pending spaces
                    mol.add_residue(resname, atoms, refresh_residues=False)
                mol.refresh_residues()

            if atom_names and all(name.strip() for name in atom_names):
                for atom in mol.atoms:
                    atom.name = atom_names[atom.id_in_mol].strip()
            else:  # don't use atom names in MAE if anyone is empty. Reassign atom names
                for residue in mol.residues:
                    for i_in_res, atom in enumerate(residue.atoms):
                        atom.name = f'{atom.symbol}{i_in_res + 1}'

            try:
                m_bond = next(block for block in f_m_ct.children if block.name == 'm_bond')
            except StopIteration:
                pass
            else:
                for i, j, order in zip(m_bond.data['i_m_from'], m_bond.data['i_m_to'], m_bond.data['i_m_order']):
                    atom1 = mol.atoms[int(i) - 1]
                    atom2 = mol.atoms[int(j) - 1]
                    mol.add_bond(atom1, atom2, int(order))

            if split_molecule == 'residue':
                pieces = mol.split_residues(consecutive=True)
            elif split_molecule == 'whole':
                pieces = [mol]
            elif split_molecule == 'atom':
                pieces = mol.split(consecutive=True)
            else:
                raise Exception(f'Invalid option for split_molecule: {split_molecule}')

            molecules.extend(pieces)

        top = Topology(molecules, cell=cell)
        top.generate_angle_dihedral_improper()

        return top

    @staticmethod
    def save_to(top, file, **kwargs):
        '''
        Save topology into an MAE file.

        Atom elements, atom names, atom types, charges, positions, formal charges, residues, bonds and unit cell will be saved.

        The MAE file will contain only one entry. All molecules in the topology will be saved in that entry.

        Parameters
        ----------
        top : Topology
        file : str
        '''
        suffix = file.split('.')[-1].lower()
        if suffix == 'mae':
            f = open(file, 'w')
        elif suffix == 'maegz':
            f = gzip.open(file, 'wt')
        else:
            raise Exception('Output must be a MAE or MAEGZ file')

        string = '''{
 s_m_m2io_version
 :::
 2.0.0
}
'''
        if top.cell.volume != 0:
            string += '''f_m_ct {
 i_m_ct_format
 r_chorus_box_ax
 r_chorus_box_ay
 r_chorus_box_az
 r_chorus_box_bx
 r_chorus_box_by
 r_chorus_box_bz
 r_chorus_box_cx
 r_chorus_box_cy
 r_chorus_box_cz
 r_pdb_PDB_CRYST1_a
 r_pdb_PDB_CRYST1_alpha
 r_pdb_PDB_CRYST1_b
 r_pdb_PDB_CRYST1_beta
 r_pdb_PDB_CRYST1_c
 r_pdb_PDB_CRYST1_gamma
 s_pdb_PDB_CRYST1_Space_Group
 :::
 2
 %.4f
 %.4f
 %.4f
 %.4f
 %.4f
 %.4f
 %.4f
 %.4f
 %.4f
 %.4f
 %.1f
 %.4f
 %.1f
 %.4f
 %.1f
 "P 1"
''' % (*flatten(top.cell.vectors), *flatten(zip(top.cell.lengths, top.cell.angles * RAD2DEG)))
        else:
            string += '''f_m_ct {
 i_m_ct_format
 :::
 2
'''
        string += ''' m_atom[%i] {
  # First column is atom index #
  i_m_atomic_number
  s_m_atom_name
  s_m_atom_type
  r_m_charge
  r_m_x_coord
  r_m_y_coord
  r_m_z_coord
  i_m_formal_charge
  i_m_residue_number
  s_m_pdb_residue_name
  :::
''' % (top.n_atom)
        for atom in top.atoms:
            string += f'  {atom.id + 1:<6d}' \
                      f' {Element(atom.symbol).number:3d}' \
                      f' {atom.name:6s}' \
                      f' {atom.type:6s}' \
                      f' {atom.charge:10.4f}' \
                      f' {atom.position[0] * 10:.4f}' \
                      f' {atom.position[1] * 10:.4f}' \
                      f' {atom.position[2] * 10:.4f}' \
                      f' {atom.formal_charge:2d}' \
                      f' {atom.residue.id + 1:4d}' \
                      f' {atom.residue.name}\n'
        string += '''  :::
 }
'''

        string += ''' m_bond[%i] {
  # First column is bond index #
  i_m_from
  i_m_to
  i_m_order
  :::
''' % (top.n_bond)
        for i, bond in enumerate(top.bonds):
            string += f'  {i + 1} {bond.atom1.id + 1} {bond.atom2.id + 1} {bond.order}\n'
        string += '''  :::
 }
}
'''

        f.write(string)
        f.close()


Topology.register_format('.mae', Mae)
Topology.register_format('.maegz', Mae)
