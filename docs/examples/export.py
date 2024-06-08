from mstk.topology import Molecule, Topology
from mstk.forcefield import ForceField, ZftTyper
from mstk.simsys import System
from mstk.wrapper import Packmol

typer = ZftTyper('alkane.zft')

butane = Molecule.from_smiles('CCCC    butane')
typer.type(butane)

# Load force field parameters from ZFF file
ff = ForceField.open('alkane.zff')

# Initialize a topology with periodic boundary condition
# For now, it contains only one butane molecule
top = Topology([butane])
top.cell.set_box([3, 3, 3])

# Assign atomic charges based on the charge parameters in the force field
ff.assign_charge(top)

# Call Packmol to build a configuration containing 100 butane molecules
packmol = Packmol(r'/path/of/packmol')
top.scale_with_packmol([100], packmol=packmol)

# Associate force field parameters to the topology
# And then export input files for simulation engines
system = System(top, ff)
system.export_gromacs()
system.export_lammps()
system.export_namd()
