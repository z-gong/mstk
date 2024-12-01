import io
import os
from .typer import Typer
from mstk.forcefield.errors import *
from mstk.topology import Bond

try:
    from rdkit.Chem import AllChem as Chem
except:
    RDKIT_IMPORTER = False
else:
    RDKIT_IMPORTER = True


class TypeDefine:
    def __init__(self, name, smarts, formal_charge=None):
        self.name = name
        self.smarts = smarts
        self.rdsmarts = None
        self.formal_charge = formal_charge  # None means don't match formal charge
        # smarts not provided means this is a fake root define
        if smarts is not None:
            try:
                self.rdsmarts = Chem.MolFromSmarts(smarts)
            except:
                raise Exception('Invalid SMARTS: %s' % smarts)

        self.children: [TypeDefine] = []
        self.parent: TypeDefine = None

    def __repr__(self):
        return '<TypeDefine: %s %s %s>' % (self.name, self.smarts, self.formal_charge)

    def add_child(self, define):
        if define not in self.children:
            self.children.append(define)
            define.parent = self


class SmartsTyper(Typer):
    '''
    SmartsTyper assign atom types with local environment defined by SMARTS and a hierarchical rule.

    A type definition file is required.

    Examples
    --------
    >>> TypeDefinition

    >>> H_1 [#1]
    >>> C_4 [#6X4]
    >>> HC  [H][CX4]
    >>> CT  [CX4;H3]
    >>> CS  [CX4;H2]
    >>> HT  [H][CX4]c1ccccc1


    >>> HierarchicalTree

    >>> H_1
    >>>     HC
    >>>         HT
    >>> C_4
    >>>     CT
    >>>     CS

    This is a sample type definition file for hydrocarbons with OPLS force field.
    Two sections are required: `TypeDefinition` and `HierarchicalTree`.
    `TypeDefinition` determines which type the atom can be by matching SMARTS.
    An atom can match several atom types defined in `TypeDefinition` section.
    e.g. Hydrogen atoms with one neighbour will match type `H_1`, and carbon atoms with four neighbour will match type `C_4`.
    Hydrogen atoms bonded with four neighboured carbon will also match type `HC`.
    If a four neighboured carbon is bonded to a benzene ring, then hydrogen atoms bonded to this carbon will also match type `HT`.
    Four neighboured carbon connected with two and three hydrogen atoms will also match `CS` and `CT`, respectively.

    After all the possible atom types are matched,
    `HierarchicalTree` determines which type the atom will finally be by performing depth first search.
    The indentation (by 4 spaces) in the `HierarchicalTree` section represents the depth of atom types.
    Therefore, if a atom matches `C_4`, but not `CT` or `CS`, it will be typed as `C_4`
    If a atom matches both `C_4` and `CT` but not `CS`, it will be typed as `CT`.
    If a atom matches both `C_4` and `CS` but not `CT`, it will be typed as `CS`.
    If a atom matches both `C_4`, `CT` and `CS` (which is impossible in this example), it will be typed as `CT`.

    For molecule propane, based on the `TypeDefinition`, all the hydrogen atoms will match atom type `H_1` and `HC`.
    All the carbon atoms will match atom type `C_4`, but the side carbon atoms will also match `CT` and the center carbon atom will also match `CS`.
    Then based on the `HierarchicalTree`, all the hydrogen atoms will be typed as `HC`.
    Side carbon atoms will be typed as `CT`, and center carbon atom will be typed as `CS`.

    The hierarchical strategy make the atom type definition extendable.

    Parameters
    ----------
    file : str or file-like object, optional
        Type definition file.

    Notes
    -----
    * SMARTS is parsed by using RDKit package. Make sure it is installed.
    * In type definition file, empty lines are ignored, and comments should start with ##.

    '''

    def __init__(self, file=None):
        if not RDKIT_IMPORTER:
            raise Exception('RDKit is required for SmartsTyper')

        self.defines: {str: TypeDefine} = {}
        self.define_root = TypeDefine('UNDEFINED', None)

        content = None
        if type(file) is str:
            if os.path.exists(file):
                with open(file) as f:
                    content = f.read()
            else:
                raise Exception(f'Typing file not found: {file}')
        elif isinstance(file, io.IOBase):
            content = file.read()
        if content is None:
            raise Exception('Argument file or string should be provided')

        self._parse(content)

    def _parse(self, content):
        lines = content.splitlines()

        section = ''
        tree_lines = []
        for line in lines:
            line = line.split('##')[0]
            if line.strip() == '':
                continue
            if line.strip() == 'TypeDefinition':
                section = 'TypeDefinition'
                continue
            if line.strip() == 'HierarchicalTree':
                section = 'HierarchicalTree'
                continue
            if section == 'TypeDefinition':
                words = line.strip().split()
                if len(words) < 2:
                    raise Exception('Atom type and SMARTS must be provided: ' + line)
                name, smarts = words[:2]
                try:
                    formal_charge = int(words[2]) if len(words) > 2 else None
                except ValueError:
                    raise Exception('Formal charge must be integer: ' + line)
                self.defines[name] = TypeDefine(name, smarts, formal_charge)
            if section == 'HierarchicalTree':
                tree_lines.append(line.rstrip())

        last_level = 0
        last_define = self.define_root
        for line in tree_lines:
            name = line.lstrip()
            indent = len(line) - len(name)
            if indent % 4 != 0:
                raise Exception('Indentation for HierarchicalTree should be 4 spaces: %s' % line)
            level = (indent) // 4 + 1
            if level == last_level + 1:
                parent = last_define
            elif level <= last_level:
                parent = last_define
                for i in range(last_level - level + 1):
                    parent = parent.parent
            else:
                raise Exception('Invalid indentation: %s' % line)
            last_level = level
            try:
                last_define = self.defines[name]
            except:
                raise Exception('Atom type not found in TypeDefinition section: %s' % line)
            parent.add_child(last_define)

    def _type_molecule(self, molecule):
        '''
        Assign atom types in the molecule with rules in type definition file.

        The :attr:`~mstk.topology.Atom.type` attribute of all atoms in the molecule will be updated.

        SmartsTyper use RDKit to do SMARTS matching, therefore the orders must be set for all the bonds in the molecule.
        Usually it means the molecule be initialized from SMILES.
        with :func:`~mstk.topology.Molecule.from_smiles` or :func:`~mstk.topology.Molecule.from_pybel`

        If an atom can not match any type by the predefined SMARTS patterns, an Exception will be raised.

        Parameters
        ----------
        molecule : Molecule
        '''
        if any(bond.order == Bond.Order.UNSPECIFIED for bond in molecule.bonds):
            raise TypingNotSupportedError(f'Typer requires all bond orders in {molecule} be set')

        possible_defines = {i: [] for i in range(molecule.n_atom)}
        for define in self.defines.values():
            for indexes in molecule.rdmol.GetSubstructMatches(define.rdsmarts, uniquify=False, maxMatches=0):
                if define.formal_charge is not None \
                        and sum(molecule.atoms[i].formal_charge for i in indexes) != define.formal_charge:
                    continue
                idx = indexes[0]
                possible_defines[idx].append(define)

        _undefined = []
        for i, atom in enumerate(molecule.atoms):
            define = self._get_deepest_define(possible_defines[i], self.define_root)
            if define is self.define_root:
                _undefined.append(atom.name)
            else:
                atom.type = define.name
        if _undefined != []:
            raise TypingUndefinedError('Definition not found for %i atoms: %s' % (
                len(_undefined), ' '.join(_undefined)))

    def _get_deepest_define(self, defines, parent: TypeDefine):
        for define in parent.children:
            if define in defines:
                return self._get_deepest_define(defines, define)
        return parent
