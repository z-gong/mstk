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


def get_omm_integrator_platform():
    integrator = mm.VerletIntegrator(0.001)
    platform = mm.Platform.getPlatformByName('Reference')
    return integrator, platform


def test_transfer_bonded_terms():
    top = Topology.open(cwd + '/../topology/files/c_3oh.msd')
    ff = ForceField.open(cwd + '/../topology/files/c_3oh.ppf')
    top.assign_charge_from_ff(ff, transfer_bci_terms=True)
    system = System(top, ff, transfer_bonded_terms=True)
    angle = next(a for a in system.topology.angles if a.name == 'C2-C3-H8')
    assert ff.get_eqt_for_angle(angle) == ('c_3', 'c_3o', 'h_1')
    aterm = system.angle_terms[id(angle)]
    assert aterm.name == 'c_3,c_3,h_1'


def test_team():
    ff = ForceField.open(cwd + '/files/10-benzene.ppf')
    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    top.assign_charge_from_ff(ff)
    system = System(top, ff)

    context = mm.Context(system.to_omm_system(), *get_omm_integrator_platform())
    context.setPositions(top.positions)

    print_energy_terms(context)

    pe = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal_mol)
    assert pytest.approx(pe, rel=0.001) == 488.5


def test_vdw_shift():
    ff = ForceField.open(cwd + '/files/10-benzene.ppf')
    ff.vdw_long_range = ff.VDW_LONGRANGE_SHIFT

    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    # top.assign_charge_from_ff(ff)
    system = System(top, ff)

    context = mm.Context(system.to_omm_system(), *get_omm_integrator_platform())
    context.setPositions(top.positions)

    print_energy_terms(context)

    pe = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal_mol)
    assert pytest.approx(pe, rel=0.001) == 494.3


def test_team_vacuum():
    ff = ForceField.open(cwd + '/files/10-benzene.ppf')
    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    top.assign_charge_from_ff(ff)
    system = System(top, ff, cell=UnitCell([0, 0, 0]))

    context = mm.Context(system.to_omm_system(), *get_omm_integrator_platform())
    context.setPositions(top.positions)

    print_energy_terms(context)

    pe = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal_mol)
    assert pytest.approx(pe, rel=0.001) == 490.8


def test_eqt_vdw():
    ff = ForceField.open(cwd + '/files/c_3ad.ppf')
    top = Topology.open(cwd + '/files/c_3ad.msd')
    top.assign_charge_from_ff(ff)
    system = System(top, ff)

    context = mm.Context(system.to_omm_system(), *get_omm_integrator_platform())
    context.setPositions(top.positions)

    print_energy_terms(context)

    pe = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal_mol)
    assert pytest.approx(pe, rel=0.001) == 46.55


def test_drude():
    ff = ForceField.open(cwd + '/../forcefield/files/CLP.ff', cwd + '/../forcefield/files/CLPol-alpha.ff')
    top = Topology.open(cwd + '/files/5-Im21-BF4-drude.lmp')
    top.generate_angle_dihedral_improper()
    # top.remove_drude_particles()
    # top.generate_drude_particles(ff)
    top.assign_charge_from_ff(ff)

    system = System(top, ff)

    context = mm.Context(system.to_omm_system(), *get_omm_integrator_platform())
    context.setPositions(top.positions)

    print_energy_terms(context)

    pe = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal_mol)
    assert pytest.approx(pe, rel=0.001) == 100.8


def test_sdk():
    ff = ForceField.open(cwd + '/../forcefield/files/SPICA_v1.zfp')
    top = Topology.open(cwd + '/files/10-SDS-20-W.lmp')
    for atom in top.atoms:
        atom.charge /= 80 ** 0.5
    system = System(top, ff)

    context = mm.Context(system.to_omm_system(), *get_omm_integrator_platform())
    context.setPositions(top.positions)

    print_energy_terms(context)

    pe = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal_mol)
    assert pytest.approx(pe, rel=0.001) == 186.0


def print_energy_terms(context):
    pe = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal_mol)
    state1 = context.getState(getEnergy=True, groups={1})
    state2 = context.getState(getEnergy=True, groups={2})
    state3 = context.getState(getEnergy=True, groups={3})
    state4 = context.getState(getEnergy=True, groups={4})
    state5 = context.getState(getEnergy=True, groups={5})
    state6 = context.getState(getEnergy=True, groups={6})
    state7 = context.getState(getEnergy=True, groups={7})
    print()
    print('E_pot      : ', pe)
    print('E_bond     : ', state1.getPotentialEnergy().value_in_unit(kcal_mol))
    print('E_angle    : ', state2.getPotentialEnergy().value_in_unit(kcal_mol))
    print('E_dihedral : ', state3.getPotentialEnergy().value_in_unit(kcal_mol))
    print('E_improper : ', state4.getPotentialEnergy().value_in_unit(kcal_mol))
    print('E_vdw      : ', state6.getPotentialEnergy().value_in_unit(kcal_mol))
    print('E_coul     : ', state5.getPotentialEnergy().value_in_unit(kcal_mol))
    print('E_drude    : ', state7.getPotentialEnergy().value_in_unit(kcal_mol))
