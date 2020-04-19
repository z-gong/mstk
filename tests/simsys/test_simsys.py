#!/usr/bin/env python3

import pytest
from mstools.topology import Topology, UnitCell
from mstools.forcefield import PpfFFSet, PaduaFFSet
from mstools.simsys import System

import os

cwd = os.path.dirname(os.path.abspath(__file__))
# top = Topology.open(cwd + '/../topology/files/benzene.zmat')
# ff = PaduaFFSet(cwd + '/../forcefield/files/oplsaa.ff')
top = Topology.open(cwd + '/files/50-Im21-BF4-drude.lmp')
ff = PaduaFFSet(cwd + '/../forcefield/files/clp.ff', cwd + '/../forcefield/files/clp-alpha.ff')

top.generate_angle_dihedral_improper()

top.remove_drude_particles()
top.assign_charge_from_forcefield(ff)

system = System(top, ff)

omm_sys = system.to_omm_system()

import simtk.openmm as mm
from simtk.openmm import app
from simtk import unit

integrator = mm.VerletIntegrator(0.001)
_platform = mm.Platform.getPlatformByName('Reference')
# _platform = mm.Platform.getPlatformByName('CPU')
# _platform = mm.Platform.getPlatformByName('OpenCL')

sim = app.Simulation(top.to_omm_topology(), omm_sys, integrator, _platform)
sim.context.setPositions(top.positions)
state0 = sim.context.getState(getEnergy=True)
state1 = sim.context.getState(getEnergy=True, groups={1})
state2 = sim.context.getState(getEnergy=True, groups={2})
state3 = sim.context.getState(getEnergy=True, groups={3})
state4 = sim.context.getState(getEnergy=True, groups={4})
state5 = sim.context.getState(getEnergy=True, groups={5})
state6 = sim.context.getState(getEnergy=True, groups={6})
state7 = sim.context.getState(getEnergy=True, groups={7})
print('E_pot:       ', state0.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
print('E_bond:      ', state1.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
print('E_angle:     ', state2.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
print('E_dihedral:  ', state3.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
print('E_improper:  ', state4.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
print('E_nb:        ', state5.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
print('E_custom_nb: ', state6.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
print('E_drude:     ', state7.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))

with open(cwd + '/files/ommsys.xml', 'w') as f:
    f.write(mm.XmlSerializer.serialize(omm_sys))
