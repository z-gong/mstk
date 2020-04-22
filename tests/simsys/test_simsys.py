#!/usr/bin/env python3

import os
import pytest
from mstools.topology import Topology, UnitCell
from mstools.forcefield import Ppf, Padua, Zfp
from mstools.simsys import System

import simtk.openmm as mm
from simtk.openmm import app
from simtk.unit import kilocalorie_per_mole as kcal_mol

cwd = os.path.dirname(os.path.abspath(__file__))


def test_team():
    ff = Ppf(cwd + '/files/10-benzene.ppf')
    top = Topology.open(cwd + '/files/10-benzene.lmp')
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
    assert pytest.approx(pe, rel=0.001) == 488.6


def test_vdw_shift():
    ff = Ppf(cwd + '/files/10-benzene.ppf')
    ff.vdw_long_range = ff.VDW_LONGRANGE_SHIFT

    top = Topology.open(cwd + '/files/10-benzene.lmp')
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


def test_drude():
    ff = Padua(cwd + '/../forcefield/files/clp.ff', cwd + '/../forcefield/files/clp-alpha.ff')
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
    assert pytest.approx(pe, rel=0.001) == 99.8


def test_sdk():
    ff = Zfp(cwd + '/../forcefield/files/SPICA_v1.zfp')
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
