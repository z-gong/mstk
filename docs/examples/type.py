from mstk.forcefield.typer import SmartsTyper
from mstk.topology import Molecule

# Construct a typing engine from type definition file
typer = SmartsTyper('alkane.smt')

# Create a molecule from SMILES string
butane = Molecule.from_smiles('CCCC')

# Assign atom types based on SMARTS pattern
typer.type(butane)
for atom in butane.atoms:
    print(atom.name, atom.type)
