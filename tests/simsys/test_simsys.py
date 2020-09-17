#!/usr/bin/env python3

import os
import sys
import pytest
from mstools.topology import Topology, UnitCell
from mstools.forcefield import ForceField
from mstools.simsys import System

import simtk.openmm as mm
from simtk.openmm import app
from simtk.unit import kilocalorie_per_mole as kcal_mol

cwd = os.path.dirname(os.path.abspath(__file__))


def test_team():
    ff = ForceField.open(cwd + '/files/10-benzene.ppf')
    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    top.assign_charge_from_ff(ff)
    system = System(top, ff)

    omm_sys = system.to_omm_system()
    integrator = mm.VerletIntegrator(0.001)
    platform = mm.Platform.getPlatformByName('Reference')
    # platform = mm.Platform.getPlatformByName('OpenCL')
    sim = app.Simulation(top.to_omm_topology(), omm_sys, integrator, platform)
    sim.context.setPositions(top.positions)

    print_energy_terms(sim)

    pe = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal_mol)
    assert pytest.approx(pe, rel=0.001) == 488.5


def test_vdw_shift():
    ff = ForceField.open(cwd + '/files/10-benzene.ppf')
    ff.vdw_long_range = ff.VDW_LONGRANGE_SHIFT

    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    # top.assign_charge_from_ff(ff)
    system = System(top, ff)

    omm_sys = system.to_omm_system()
    integrator = mm.VerletIntegrator(0.001)
    platform = mm.Platform.getPlatformByName('Reference')
    # platform = mm.Platform.getPlatformByName('OpenCL')
    sim = app.Simulation(top.to_omm_topology(), omm_sys, integrator, platform)
    sim.context.setPositions(top.positions)

    print_energy_terms(sim)

    pe = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal_mol)
    assert pytest.approx(pe, rel=0.001) == 494.3


def test_team_vacuum():
    ff = ForceField.open(cwd + '/files/10-benzene.ppf')
    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    top.assign_charge_from_ff(ff)
    system = System(top, ff, cell=UnitCell([0, 0, 0]))

    omm_sys = system.to_omm_system()
    integrator = mm.VerletIntegrator(0.001)
    platform = mm.Platform.getPlatformByName('Reference')
    sim = app.Simulation(top.to_omm_topology(), omm_sys, integrator, platform)
    sim.context.setPositions(top.positions)

    print_energy_terms(sim)

    pe = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal_mol)
    assert pytest.approx(pe, rel=0.001) == 490.8


def test_drude():
    ff = ForceField.open(cwd + '/../forcefield/files/CLP.ff', cwd + '/../forcefield/files/CLPol-alpha.ff')
    top = Topology.open(cwd + '/files/5-Im21-BF4-drude.lmp')
    top.generate_angle_dihedral_improper()
    # top.remove_drude_particles()
    # top.generate_drude_particles(ff)
    top.assign_charge_from_ff(ff)

    system = System(top, ff)

    omm_sys = system.to_omm_system()
    integrator = mm.VerletIntegrator(0.001)
    platform = mm.Platform.getPlatformByName('Reference')
    # platform = mm.Platform.getPlatformByName('OpenCL')
    sim = app.Simulation(top.to_omm_topology(), omm_sys, integrator, platform)
    sim.context.setPositions(top.positions)

    print_energy_terms(sim)

    pe = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal_mol)
    assert pytest.approx(pe, rel=0.001) == 100.8


def test_sdk():
    ff = ForceField.open(cwd + '/../forcefield/files/SPICA_v1.zfp')
    top = Topology.open(cwd + '/files/10-SDS-20-W.lmp')
    for atom in top.atoms:
        atom.charge /= 80 ** 0.5

    system = System(top, ff)
    omm_top = top.to_omm_topology()
    omm_sys = system.to_omm_system()

    integrator = mm.VerletIntegrator(0.001)

    platform = mm.Platform.getPlatformByName('Reference')
    # platform = mm.Platform.getPlatformByName('OpenCL')
    sim = app.Simulation(omm_top, omm_sys, integrator, platform)
    sim.context.setPositions(top.positions)

    print_energy_terms(sim)

    pe = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal_mol)
    assert pytest.approx(pe, rel=0.001) == 186.0


def print_energy_terms(sim):
    pe = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal_mol)
    state1 = sim.context.getState(getEnergy=True, groups={1})
    state2 = sim.context.getState(getEnergy=True, groups={2})
    state3 = sim.context.getState(getEnergy=True, groups={3})
    state4 = sim.context.getState(getEnergy=True, groups={4})
    state5 = sim.context.getState(getEnergy=True, groups={5})
    state6 = sim.context.getState(getEnergy=True, groups={6})
    state7 = sim.context.getState(getEnergy=True, groups={7})
    print()
    print('E_pot      : ', pe)
    print('E_bond     : ', state1.getPotentialEnergy().value_in_unit(kcal_mol))
    print('E_angle    : ', state2.getPotentialEnergy().value_in_unit(kcal_mol))
    print('E_dihedral : ', state3.getPotentialEnergy().value_in_unit(kcal_mol))
    print('E_improper : ', state4.getPotentialEnergy().value_in_unit(kcal_mol))
    print('E_vdw      : ', state6.getPotentialEnergy().value_in_unit(kcal_mol))
    print('E_coul     : ', state5.getPotentialEnergy().value_in_unit(kcal_mol))
    print('E_drude    : ', state7.getPotentialEnergy().value_in_unit(kcal_mol))
