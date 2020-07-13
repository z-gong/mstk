from mstools.topology import Molecule, Topology, UnitCell
from mstools.forcefield import ForceField, ZftTyper
from mstools.simsys import System
from mstools.wrapper import Packmol

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

ff = ForceField.open('alkane.zfp')

top = Topology([butane], cell=UnitCell([3, 3, 3]))
top.assign_charge_from_ff(ff)
packmol = Packmol(r'/Users/zgong/miniconda3/bin/packmol')
top.scale_with_packmol([100], packmol=packmol)

system = System(top, ff)
system.export_gromacs()
system.export_lammps()
