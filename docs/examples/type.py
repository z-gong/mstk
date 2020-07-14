from mstools.forcefield import ZftTyper
from mstools.topology import Molecule

# Construct a typing engine from type definition
definition = '''
TypeDefinition

h_1    [H][CX4]         0
c_4    [CX4]            0
c_4h2  [CX4;H2]         0
c_4h3  [CX4;H3]         0

HierarchicalTree

h_1
c_4
    c_4h2
    c_4h3
'''
typer = ZftTyper(content=definition)

# Initialize a molecule from SMILES string
butane = Molecule.from_smiles('CCCC    butane')

# Assign atom types based on SMARTS pattern
typer.type_molecule(butane)
for atom in butane.atoms:
    print(atom.name, atom.type)
