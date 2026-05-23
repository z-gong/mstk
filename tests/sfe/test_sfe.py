"""
Uses the existing 10-benzene system to test the alchemical SFE module (mstk.sfe).
"""

import os
import openmm.openmm as mm
import pytest
from openmm import app, unit
from mstk.topology import Topology
from mstk.forcefield import ForceField
from mstk.simsys import System
from mstk.sfe import SFEManager

cwd = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(cwd, '..', 'simsys', 'files')


@pytest.fixture
def benzene_system():
    top = Topology.open(os.path.join(DATA_DIR, '10-benzene.lmp'), improper_center=3)
    ff = ForceField.open(os.path.join(DATA_DIR, '10-benzene.ppf'))
    ff.assign_charge(top)
    return System(top, ff)


@pytest.fixture
def alchemical_context(benzene_system):
    sfe = SFEManager(benzene_system, 0)
    platform = mm.Platform.getPlatformByName('Reference')
    integrator = mm.VerletIntegrator(0.001)
    omm_top = benzene_system.topology.to_omm_topology()
    sim = app.Simulation(omm_top, sfe.omm_system, integrator, platform)
    sim.context.setPositions(benzene_system.topology.positions)
    return sim.context, sfe


class TestLambdaSchedule:
    def test_endpoints(self):
        schedule = SFEManager.default_lambda_schedule(16)
        assert schedule[0] == (1.0, 1.0)
        assert schedule[-1] == (0.0, 0.0)

    def test_coul_first(self):
        schedule = SFEManager.default_lambda_schedule(16)
        for lam_coul, lam_vdw in schedule:
            if lam_coul > 0:
                assert lam_vdw == 1.0

    def test_vdw_after_coul(self):
        schedule = SFEManager.default_lambda_schedule(16)
        for lam_coul, lam_vdw in schedule:
            if lam_vdw < 1.0:
                assert lam_coul == 0.0

    def test_monotonic(self):
        schedule = SFEManager.default_lambda_schedule(16)
        coul_values = [lc for lc, lv in schedule if lv == 1.0]
        vdw_values = [lv for lc, lv in schedule if lc == 0.0]
        assert coul_values == sorted(coul_values, reverse=True)
        assert vdw_values == sorted(vdw_values, reverse=True)

    def test_n_windows(self):
        for n in [8, 12, 16, 20, 24]:
            schedule = SFEManager.default_lambda_schedule(n)
            assert len(schedule) == n


class TestCoupledEndState:
    def test_energy_offset_small(self, benzene_system, alchemical_context):
        """At lambda=(1,1), alchemical system differs from reference by a small LRC offset."""
        context, sfe = alchemical_context
        sfe.set_lambda_state(context, 0)
        pe_alch = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoules_per_mole)

        omm_ref = benzene_system.to_omm_system()
        platform = mm.Platform.getPlatformByName('Reference')
        sim_ref = app.Simulation(
            benzene_system.topology.to_omm_topology(), omm_ref,
            mm.VerletIntegrator(0.001), platform)
        sim_ref.context.setPositions(benzene_system.topology.positions * unit.nanometer)
        pe_ref = sim_ref.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoules_per_mole)

        offset = pe_alch - pe_ref
        assert abs(offset) < 1.0

    def test_bonded_energies_unchanged(self, benzene_system, alchemical_context):
        """Bonded energies must be identical to reference."""
        context, sfe = alchemical_context
        sfe.set_lambda_state(context, 0)

        omm_ref = benzene_system.to_omm_system()
        platform = mm.Platform.getPlatformByName('Reference')
        sim_ref = app.Simulation(
            benzene_system.topology.to_omm_topology(), omm_ref,
            mm.VerletIntegrator(0.001), platform)
        sim_ref.context.setPositions(benzene_system.topology.positions * unit.nanometer)

        for group in [1, 2, 3, 4]:
            pe_alch = context.getState(getEnergy=True, groups={group}).getPotentialEnergy()
            pe_ref = sim_ref.context.getState(getEnergy=True, groups={group}).getPotentialEnergy()
            assert abs(pe_alch - pe_ref) < 1e-6 * unit.kilojoules_per_mole


class TestDecoupledEndState:
    def test_alchemical_lj_zero(self, alchemical_context):
        """At lambda_vdw=0, alchemical LJ force (group 11) should be zero."""
        context, sfe = alchemical_context
        sfe.set_lambda_state(context, 15)
        pe_alch_lj = context.getState(getEnergy=True, groups={11}).getPotentialEnergy()
        assert abs(pe_alch_lj.value_in_unit(unit.kilojoules_per_mole)) < 1e-6

    def test_energy_matches_separated(self, benzene_system, alchemical_context):
        """At lambda=(0,0), energy should equal E_env + E_M (up to LRC offset)."""
        context, sfe = alchemical_context
        sfe.set_lambda_state(context, 15)
        pe_decoupled = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoules_per_mole)

        top = benzene_system.topology
        ff = benzene_system.ff
        platform = mm.Platform.getPlatformByName('Reference')

        top_env = Topology(top.molecules[1:])
        top_env.cell = top.cell
        system_env = System(top_env, ff)
        omm_env = system_env.to_omm_system()
        sim_env = app.Simulation(
            top_env.to_omm_topology(), omm_env,
            mm.VerletIntegrator(0.001), platform)
        sim_env.context.setPositions(top_env.positions)
        pe_env = sim_env.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoules_per_mole)

        top_m = Topology([top.molecules[0]])
        top_m.cell = top.cell
        system_m = System(top_m, ff)
        omm_m = system_m.to_omm_system()
        sim_m = app.Simulation(
            top_m.to_omm_topology(), omm_m,
            mm.VerletIntegrator(0.001), platform)
        sim_m.context.setPositions(top_m.positions)
        pe_m = sim_m.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoules_per_mole)

        offset = pe_decoupled - (pe_env + pe_m)
        assert abs(offset) < 1.0


class TestEnergySmoothing:
    def test_no_discontinuity(self, alchemical_context):
        """No large jumps between consecutive windows."""
        context, sfe = alchemical_context

        energies = []
        for i in range(16):
            sfe.set_lambda_state(context, i)
            pe = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
                unit.kilojoules_per_mole)
            energies.append(pe)

        for i in range(len(energies) - 1):
            dU = abs(energies[i + 1] - energies[i])
            assert dU < 50.0, f"Discontinuity at window {i}->{i+1}: dU={dU:.2f} kJ/mol"

    def test_du_self_is_zero(self, alchemical_context):
        """dU to the same state must be exactly zero."""
        context, sfe = alchemical_context

        sfe.set_lambda_state(context, 5)
        U1 = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoules_per_mole)
        sfe.set_lambda_state(context, 5)
        U2 = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoules_per_mole)
        assert abs(U2 - U1) < 1e-10


class TestGuards:
    def test_reject_vdw_shift(self, benzene_system):
        benzene_system.ff.vdw_long_range = ForceField.VDW_LONGRANGE_SHIFT
        with pytest.raises(ValueError, match='shift'):
            SFEManager(benzene_system, 0)

    def test_reject_polarizable(self, benzene_system):
        benzene_system.ff.polar_terms['fake'] = None
        with pytest.raises(ValueError, match='[Pp]olarizable'):
            SFEManager(benzene_system, 0)

    def test_reject_virtual_site(self, benzene_system):
        benzene_system.ff.virtual_site_terms['fake'] = None
        with pytest.raises(ValueError, match='[Vv]irtual site'):
            SFEManager(benzene_system, 0)
