from mstools.forcefield import ZftTyper
from mstools.topology import Molecule, Topology

definition = '''
TypeDefinition

H_1 [#1]              0
C_3 [#6X3]            0
C_4 [#6X4]            0
HC  [H][CX4]         0
CT  [CX4;H3]         0
CS  [CX4;H2]         0
HA  [H]c1ccccc1      0
CA  c1ccccc1         0


HierarchicalTree

H_1
    HC
    HA
C_3
    CA
C_4
    CT
    CS
'''

typer = ZftTyper(content=definition)

butane = Molecule.from_smiles('CCCC    butane')
typer.type_molecule(butane)
for atom in butane.atoms:
    print(atom.name, atom.type)

benzene = Molecule.from_smiles('c1ccccc1 benzene')
typer.type_molecule(benzene)
for atom in benzene.atoms:
    print(atom.name, atom.type)
