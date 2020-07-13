from mstools.forcefield import ZftTyper
from mstools.topology import Molecule

definition = '''
TypeDefinition

h_1    [H][CX4]         0
c_4h3  [CX4;H3]         0
c_4h2  [CX4;H2]         0

HierarchicalTree

h_1
c_4h3
c_4h2
'''

typer = ZftTyper(content=definition)

butane = Molecule.from_smiles('CCCC    butane')
typer.type_molecule(butane)
for atom in butane.atoms:
    print(atom.name, atom.type)
