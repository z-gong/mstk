# mstk

_A toolkit to make molecular simulation less painful_

`mstk` is a Python toolkit designed to streamline the setup and management of molecular simulations, particularly geared
towards non-biological applications in material science and chemical engineering. It simplifies atom typing, force field
parameter assignment, and input file generation for major simulation engines.

## Features

* Assign atom types and force field parameters based on local chemical environment
* Generate input files for LAMMPS, GROMACS, NAMD and OpenMM
* Set up solvation free energy simulation
* Support Drude polarizable model, coarse-grained model, virtual site, linear angle...
* Read/write common topology (PSF, ZMAT, PDB, LAMMPS, ...) and trajectory (GRO, XTC, DCD, LAMMPS, ...) files
* Access local and remote job schedulers like Slurm

## Installation

```
conda install -c conda-forge numpy pandas rdkit openmm chemfiles packmol pymbar-core
pip install mstk
```

## Quick Example

**Example 1** Build a liquid mixture of 100 benzene and 1000 water molecules, ready for simulation.

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


**Example 2** Calculate the hydration free energy of benzene in water

- 2.1 Build a system containing 1 benzene in 1000 water, following the above example

- 2.2 Use the `mstk sfe run` command to perform alchemical simulation under different lambda windows
```
mstk sfe run -p _namd-top.psf -c _gmx-conf.gro -f primitive.zff --molid 0 -i %WINDOW% --nwindow 16
```
Run at 16 different windows by change %WINDOW% from 0 to 15.
A series of dU values will be saved in dU_%WINDOW%.csv

- 2.3 Use the `mstk sfe mbar` command to analyse free energy change
```
mstk sfe mbar *.csv
```

## Documentation

https://mstk.readthedocs.io/en/latest/index.html

## TODO

- [ ] Refactor algorithms across `topology` and `analyzer` modules
- [ ] Separate `ommhelper` module into its own package
- [ ] Remove `analyzer` module
- [ ] Revise the implementation of Drude polarization and virtual sites
