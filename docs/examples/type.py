from mstk.forcefield.typer import ZftTyper
from mstk.topology import Molecule

# Construct a typing engine from type definition file
typer = ZftTyper('alkane.zft')

# Create a molecule from SMILES string
butane = Molecule.from_smiles('CCCC')

# Assign atom types based on SMARTS pattern
typer.type(butane)
for atom in butane.atoms:
    print(atom.name, atom.type)
