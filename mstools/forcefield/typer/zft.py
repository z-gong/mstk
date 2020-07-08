from .typer import Typer
from ..errors import *

try:
    import openbabel as ob
    import pybel
except:
    OPENBABEL_IMPORTED = False
else:
    OPENBABEL_IMPORTED = True


class TypeDefine():
    def __init__(self, name, smarts, charge):
        self.name = name
        self.obsmarts = None
        self.charge = charge
        if smarts is not None:
            self.obsmarts = ob.OBSmartsPattern()
            if not self.obsmarts.Init(smarts):
                raise Exception('Invalid SMARTS: %s' % smarts)

        self.children: [TypeDefine] = []
        self.parent: TypeDefine = None

    def __repr__(self):
        return '<TypeDefine: %s>' % self.name

    def add_child(self, define):
        if define not in self.children:
            self.children.append(define)
            define.parent = self


class ZftTyper(Typer):
    '''
    ZftTyper assign atom types with local environment defined by SMARTS and a hierarchical rule.

    A type definition file is required.

    Examples
    --------
    >>> TypeDefinition

    >>> H_1 [H]              0
    >>> C_4 [#6X4]           0
    >>> HC  [H][CX4]         0
    >>> CT  [CX4;H3]         0
    >>> CS  [CX4;H2]         0
    >>> HT  [H][CX4]c1ccccc1 0


    >>> HierarchicalTree

    >>> H_1
    >>>     HC
    >>>         HT
    >>> C_4
    >>>     CT
    >>>     CS

    This is a sample type definition file for hydrocarbons with OPLS force field.
    Two sections are required: `TypeDefinition` and `HierarchicalTree`.
    `TypeDefinition` determines which type the atom can be by matching SMARTS and formal charge.
    Currently, charge information is ignored, but it should be provided in the file.
    A atom can match several different atom types defined in `TypeDefinition` section.
    e.g. Hydrogen atoms with one neighbour will match type `H_1`, and carbon atoms with four neighbour will match type `C_4`.
    Hydrogen atoms bonded with four neighboured carbon will also match type `H_C`.
    Four neighboured carbon connected with three hydrogen atoms will also match `CT`.

    After all the possible atom types are matched,
    `HierarchicalTree` determines which type the atom will finally be by performing depth first search.
    The indentation (by 4 spaces) in the `HierarchicalTree` section represents the depth of atom types.
    Therefore, if a atom matches `C_4`, but not `CT` or `CS`, it will be typed as `C_4`
    If a atom matches both `C_4` and `CT` but not `CS`, it will be typed as `CT`.
    If a atom matches both `C_4` and `CS` but not `CT`, it will be typed as `CS`.
    If a atom matches both `C_4`, `CT` and `CS` (which is impossible in this example), it will be typed as `CT`.

    For molecule propane, based on the `TypeDefinition`, all the hydrogen atoms will match atom type `H_1` and `HC`.
    All of the carbon atoms will match atom type `C_4`, but the side carbon atoms will also match `CT` and the center carbon atom will also match `CS`.
    Then based on the `HierarchicalTree`, all of the hydrogen atoms will be typed as `HC`.
    Side carbon atoms will be typed as `CT`, and center carbon atom will be typed as `CS`.

    The hierarchical strategy make the atom type definition extendable.

    Parameters
    ----------
    file : str
        Name of type definition file

    Notes
    -----
    * SMARTS is parsed by using OpenBabel package. Therefore `pybel` module should be installed.
    * In type definition file, empty lines are ignored, and comments should start with ##.

    '''

    def __init__(self, file):
        if not OPENBABEL_IMPORTED:
            raise Exception('Cannot import openbabel')

        self.defines: {str: TypeDefine} = {}
        self.define_root = TypeDefine('UNDEFINED', None, 0)
        self._file = file
        self._parse(file)

    def __repr__(self):
        return '<ZftTyper: %s>' % self._file

    def _parse(self, file):
        with open(file)  as f:
            lines = f.read().splitlines()

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
                try:
                    name, smarts, charge = line.strip().split()
                except:
                    raise Exception('smarts and charge should be provided: %s' % line)
                self.defines[name] = TypeDefine(name, smarts, charge)
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

    def type_molecule(self, molecule):
        '''
        Assign atom types in the molecule with rules in type definition file.

        The :attr:`~mstools.topology.Atom.type` attribute of all atoms in the molecule will be updated.

        ZftTyper use OpenBabel to do SMARTS matching,
        therefore it expects :attr:`~mstools.topology.Molecule.obmol` attribute to be available in the molecule.
        Usually it means the molecule should be initialized from SMILES or pybel Molecule
        with :func:`~mstools.topology.Molecule.from_smiles` or :func:`~mstools.topology.Molecule.from_pybel`

        If :attr:`~mstools.topology.Molecule.obmol` attribute is None, an Exception will be raised.
        If an atom can not match any type by the predefined SMARTS patterns, an Exception will be raised.

        Parameters
        ----------
        molecule : Molecule
        '''
        if molecule.obmol is None:
            raise TypingNotSupportedError('obmol attribute not found for %s' % str(molecule))

        possible_defines = {i: [] for i in range(molecule.n_atom)}
        for define in self.defines.values():
            obsmarts = define.obsmarts
            obsmarts.Match(molecule.obmol)
            results = list(obsmarts.GetMapList())
            for indexes in results:
                idx = indexes[0] - 1
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
