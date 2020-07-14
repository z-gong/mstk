from mstools.topology import Molecule, Topology, UnitCell
from mstools.forcefield import ForceField, ZftTyper
from mstools.simsys import System
from mstools.wrapper import Packmol

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

butane = Molecule.from_smiles('CCCC    butane')
typer.type_molecule(butane)

# Load force field parameters from ZFP file
ff = ForceField.open('alkane.zfp')

# Initialize a topology with periodic boundary condition
# For now, it contains only one butane molecule
top = Topology([butane], cell=UnitCell([3, 3, 3]))

# Assign atomic charges based on the charge parameters in the force field
top.assign_charge_from_ff(ff)

# Call Packmol to build a configuration containing 100 butane molecules
packmol = Packmol(r'/path/of/packmol')
top.scale_with_packmol([100], packmol=packmol)

# Associate force field parameters to the topology
# And then export input files for simulation engines
system = System(top, ff)
system.export_gromacs()
system.export_lammps()
system.export_charmm()
