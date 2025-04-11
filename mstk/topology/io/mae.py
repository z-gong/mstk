import gzip
import shlex
from mstk.topology import Atom, Bond, Molecule, Topology, UnitCell
from mstk.chem.element import Element
from mstk.chem.constant import DEG2RAD


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
    Generate Topology from MAE or MAEGZ file of Schrödinger.

    The atoms, positions, bonds are parsed.
    Formal charges, residues will be parsed if avaiable.
    The first available unit cell from the structure entryies will be parsed.
    The atom names are ignored, because atom name in MAE can be empty.
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

            mol = Molecule(name=f_m_ct.data['s_m_title'].replace(' ', '_'))
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

            # don't use atom names in MAE. They could be empty. Reassign atom names
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


Topology.register_format('.mae', Mae)
Topology.register_format('.maegz', Mae)
