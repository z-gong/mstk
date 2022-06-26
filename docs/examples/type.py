import io
from mstk.forcefield.typer import ZftTyper
from mstk.topology import Molecule

# Construct a typing engine from type definition
definition = '''
TypeDefinition

h_1    [H][CX4]
c_4    [CX4]
c_4h2  [CX4;H2]
c_4h3  [CX4;H3]

HierarchicalTree

h_1
c_4
    c_4h2
    c_4h3
'''
typer = ZftTyper(io.StringIO(definition))

# Create a molecule from SMILES string
butane = Molecule.from_smiles('CCCC')

# Assign atom types based on SMARTS pattern
typer.type(butane)
for atom in butane.atoms:
    print(atom.name, atom.type)
