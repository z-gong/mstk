import json
from .typer import Typer

try:
    import openbabel as ob
    import pybel
except:
    OPENBABEL_IMPORTED = True
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
    def __init__(self, file):
        super().__init__()
        self.defines: {str: TypeDefine} = {}
        self.define_root = TypeDefine('Root', None, 0)
        self._file = file
        self._parse(file)

    def __repr__(self):
        return '<ZftTyper: %s>' % self._file

    def _parse(self, file):
        if not OPENBABEL_IMPORTED:
            raise Exception('Cannot import openbabel')

        with open(file)  as f:
            lines = f.read().splitlines()

        string_tree = ''
        section = ''
        for line in lines:
            line = line.split('##')[0].strip()
            if line == '':
                continue
            if line == 'TypeDefinition':
                section = line
                continue
            if line == 'HierarchicalTree':
                section = line
                continue
            if section == 'TypeDefinition':
                name, smarts, charge = line.strip().split()
                self.defines[name] = TypeDefine(name, smarts, charge)
            if section == 'HierarchicalTree':
                string_tree += line

        try:
            tree = json.loads(string_tree)
        except Exception as e:
            raise Exception('Cannot parse HierarchicalTree, make sure it\'s valid json: %s' % e)

        self._parse_tree(self.define_root, tree)

    def _parse_tree(self, parent: TypeDefine, d: {}):
        for k, v in d.items():
            define = self.defines[k]
            parent.add_child(define)
            if v is None:
                continue
            else:
                self._parse_tree(define, v)

    def type_molecule(self, molecule):
        obmol = molecule._obmol
        if obmol is None:
            raise Exception('Cannot type this molecule: obmol not found')

        possible_defines = {i: [] for i in range(molecule.n_atom)}
        for define in self.defines.values():
            obsmarts = define.obsmarts
            obsmarts.Match(obmol)
            results = list(obsmarts.GetMapList())
            for indexes in results:
                idx = indexes[0] - 1
                possible_defines[idx].append(define)

        for k, v in possible_defines.items():
            print(molecule.atoms[k], v)

        for i in range(molecule.n_atom):
            define = self._get_deepest_define(possible_defines[i], self.define_root)
            molecule.atoms[i].type = define.name

    def _get_deepest_define(self, defines, parent: TypeDefine):
        for define in parent.children:
            if define in defines:
                return self._get_deepest_define(defines, define)
        return parent
