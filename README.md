# mstk

_A toolkit to make molecular simulation less painful_

`mstk` is a Python toolkit designed to streamline the setup and management of molecular simulations, particularly geared
towards non-biological applications in material science and chemical engineering. It simplifies atom typing, force field
parameter assignment, and input file generation for major simulation engines.

## Features

* Assign atom types and force field parameters based on local chemical environment
* Generate input files for LAMMPS, GROMACS, NAMD and OpenMM
* Support Drude polarizable model, coarse-grained model, virtual site, linear angle...
* Read/write common topology (PSF, ZMAT, PDB, LAMMPS, ...) and trajectory (GRO, XTC, DCD, LAMMPS, ...) files
* Access local and remote job schedulers like Slurm

## Installation

```
conda install -c conda-forge numpy pandas rdkit openmm chemfiles packmol
pip install mstk
```

## Quick Example

Build a liquid mixture of 100 benzene and 1000 water molecules, ready for simulation.

```python
from mstk.topology import Molecule, Topology
from mstk.forcefield import ForceField, Typer
from mstk.simsys import System
from mstk.wrapper import Packmol

# Create topology from SMILES
benzene = Molecule.from_smiles('c1ccccc1')
water = Molecule.from_smiles('O')
top = Topology([benzene, water])

# Assign atom types as defined in `data/forcefield/primitive.smt`
typer = Typer.open('primitive.smt')
typer.type(top)

# Build a bulk liquid simulation box with Packmol
packmol = Packmol('packmol')
top.cell.set_box([4.0, 4.0, 4.0])
top.scale_with_packmol([100, 1000], packmol=packmol)

# Assign force field parameters as defined in `data/forcefield/primitive.zff`
ff = ForceField.open('primitive.zff')
ff.assign_charge(top)
system = System(top, ff)

# Generate input files for LAMMPS, GROMACS and NAMD
system.export_lammps()
system.export_gromacs()
system.export_namd()

# Generate OpenMM system and topology
omm_sys = system.to_omm_system()
omm_top = top.to_omm_topology()
```

__Important Note__:
*The __primitive__ atom typing and force field used above are for demonstration purpose only and are __not__ well
optimized for production use.*
*For reliable simulation, please prepare you own `smt` and `zff` files or get them from a validated source.*

## Documentation

https://mstk.readthedocs.io/en/latest/index.html

## TODO

- [ ] Refactor algorithms across `topology` and `analyzer` modules
- [ ] Separate `ommhelper` module into its own package
- [ ] Remove `analyzer` module
- [ ] Revise the implementation of Drude polarization and virtual sites
