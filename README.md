# mstk

A toolkit to make molecular simulation less painful

## Capabilities

* Assign atom types and force field parameters based on local chemical environment
* Generate input files for simulation engines like LAMMPS, GROMACS, NAMD, OpenMM
* Support Drude polarizable model, coarse-grained model, virtual site, linear angle...
* Read/write common topology files (LAMMPS, PSF, PDB, ZMAT etc...)
* Read/write common trajectory files (LAMMPS, GRO, DCD, XTC, etc...)
* Access local and remote job schedulers like Slurm, Torque

## Examples

The following code builds a liquid mixture of 100 benzene and 1000 water molecules, ready for simulation.

```python
from mstk.topology import Molecule, Topology
from mstk.forcefield import ForceField, ZftTyper
from mstk.simsys import System
from mstk.wrapper import Packmol

# Create molecule structures from SMILES
benzene = Molecule.from_smiles('c1ccccc1')
water = Molecule.from_smiles('O')
top = Topology([benzene, water])

# Assign atom type as defined in `data/forcefield/primitive.zft`
typer = ZftTyper('primitive.zft')
typer.type(top)

# Build a bulk liquid simulation box with Packmol
packmol = Packmol('/path/to/packmol')
top.cell.set_box([4.0, 4.0, 4.0])
top.scale_with_packmol([100, 1000], packmol=packmol)

# Assign force field as defined in `data/forcefield/primitive.zff`
ff = ForceField.open('primitive.zff')
ff.assign_charge(top)
system = System(top, ff)

# generate input files for LAMMPS, GROMACS and NAMD
system.export_lammps()
system.export_gromacs()
system.export_namd()

# generate OpenMM system and topology
omm_sys = system.to_omm_system()
omm_top = top.to_omm_topology()
```

## Documentation

https://mstk.readthedocs.io/en/latest/index.html

## TODO

- [ ] Design a file format to store all the information in a topology
- [ ] Take bond order into consideration for force field assignment
- [ ] Re-organize algorithms scattered in topology and analyzer modules
- [ ] Remove dependency on `mdtraj`

## Known issue

`mstk` use `chemfiles` to write `XTC` trajectory. However, lastest `chemfiles` fails writing binary trajectory format
under WSL. If you are working under WSL, please install `chemfiles 0.10.2`

```
conda install chemfiles-lib=0.10.2 chemfiles-python=0.10.2
```
