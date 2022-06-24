#!/usr/bin/env python3

import os
import itertools
import tempfile
import filecmp
import pytest
from openmm import openmm as mm
from openmm import unit
from mstools.topology import Topology, Molecule
from mstools.analyzer.energy_kernels import *

cwd = os.path.dirname(os.path.abspath(__file__))

top = Topology.open(cwd + '/../topology/files/CH3NH2.pdb')


def calc_energy_with_omm(force, positions):
    system = mm.System()
    for i in range(len(positions)):
        system.addParticle(1)
    system.addForce(force)

    integrator = mm.VerletIntegrator(0.001)
    platform = mm.Platform.getPlatformByName('Reference')
    context = mm.Context(system, integrator, platform)
    context.setPositions(positions)

    state = context.getState(getEnergy=True, getForces=True)
    return state.getPotentialEnergy().value_in_unit_system(unit.md_unit_system), \
           state.getForces(asNumpy=True).value_in_unit_system(unit.md_unit_system)


def test_harmonic_bond_kernel():
    kernel = HarmonicBondKernel(top.positions,
                                [bond.id_atoms for bond in top.bonds],
                                [[0.1 + 0.01 * i, 50 + 5 * i] for i in range(top.n_bond)],
                                )
    r, energy, forces = kernel.evaluate()
    print()
    print(r)
    print(energy)
    print(sum(energy))
    print(forces)

    force = mm.HarmonicBondForce()
    force.setForceGroup(1)
    for i, bond in enumerate(top.bonds):
        force.addBond(*bond.id_atoms, 0.1 + 0.01 * i, 100 + 10 * i)
    e, f = calc_energy_with_omm(force, top.positions)

    assert pytest.approx(sum(energy), rel=1E-6) == e
    assert pytest.approx(forces, rel=1E-6) == f


def test_harmonic_angle_kernel():
    kernel = HarmonicAngleKernel(top.positions,
                                 [angle.id_atoms for angle in top.angles],
                                 [[1 + 0.1 * i, 50 + 5 * i] for i in range(top.n_angle)],
                                 )
    theta, energy, forces = kernel.evaluate()
    print()
    print(theta * 180 / np.pi)
    print(energy)
    print(sum(energy))
    print(forces)

    force = mm.HarmonicAngleForce()
    force.setForceGroup(2)
    for i, angle in enumerate(top.angles):
        force.addAngle(*angle.id_atoms, 1 + 0.1 * i, 100 + 10 * i)
    e, f = calc_energy_with_omm(force, top.positions)

    assert pytest.approx(sum(energy), rel=1E-6) == e
    assert pytest.approx(forces, rel=1E-6) == f


def test_opls_torsion_kernel():
    kernel = OplsTorsionKernel(top.positions,
                               [dihedral.id_atoms for dihedral in top.dihedrals],
                               [[1 + 0.1 * i, 2 + 0.1 * i, 3 + 0.1 * i, 4 + 0.1 * i]
                                for i in range(top.n_dihedral)],
                               )
    phi, energy, forces = kernel.evaluate()
    print()
    for dihedral, val in zip(top.dihedrals, phi):
        print([i + 1 for i in dihedral.id_atoms], val * 180 / np.pi)
    print(energy)
    print(sum(energy))
    print(forces)

    force = mm.PeriodicTorsionForce()
    force.setForceGroup(3)
    for i, dihedral in enumerate(top.dihedrals):
        for k in range(1, 5):
            force.addTorsion(*dihedral.id_atoms, k, ((k - 1) % 2) * np.pi, k + 0.1 * i)
    e, f = calc_energy_with_omm(force, top.positions)

    assert pytest.approx(sum(energy), rel=1E-6) == e
    assert pytest.approx(forces, rel=1E-6) == f


def test_harmonic_torsion_kernel():
    kernel = HarmonicTorsionKernel(top.positions,
                                   [dihedral.id_atoms for dihedral in top.dihedrals],
                                   [[1 + 0.1 * i, 1 + 0.05 * i] for i in range(top.n_dihedral)]
                                   )
    phi, energy, forces = kernel.evaluate()
    print()
    print(phi * 180 / np.pi)
    print(energy)
    print(sum(energy))
    print(forces)

    force = mm.CustomTorsionForce('0.5*k*min(2*pi-dtheta, dtheta)^2;'
                                  'dtheta=abs(theta-theta0);'
                                  'pi=3.1415926535')
    force.addPerTorsionParameter('theta0')
    force.addPerTorsionParameter('k')
    force.setForceGroup(4)
    for i, dihedral in enumerate(top.dihedrals):
        force.addTorsion(*dihedral.id_atoms, [1 + 0.1 * i, 2 + 0.1 * i])

    e, f = calc_energy_with_omm(force, top.positions)
    print(e)
    print(f)

    assert pytest.approx(sum(energy), rel=1E-6) == e
    assert pytest.approx(forces, rel=1E-6) == f


def test_constrained_torsion_kernel():
    kernel = ConstrainedTorsionKernel(top.positions,
                                      [dihedral.id_atoms for dihedral in top.dihedrals],
                                      [[1 + 0.1 * i, 2 + 0.1 * i] for i in range(top.n_dihedral)]
                                      )
    phi, energy, forces = kernel.evaluate()
    print()
    print(phi * 180 / np.pi)
    print(energy)
    print(sum(energy))
    print(forces)

    force = mm.CustomTorsionForce('k*(1-cos(theta-theta0))')
    force.addPerTorsionParameter('theta0')
    force.addPerTorsionParameter('k')
    force.setForceGroup(4)
    for i, dihedral in enumerate(top.dihedrals):
        force.addTorsion(*dihedral.id_atoms, [1 + 0.1 * i, 2 + 0.1 * i])

    e, f = calc_energy_with_omm(force, top.positions)
    print(e)
    print(f)

    assert pytest.approx(sum(energy), rel=1E-6) == e
    assert pytest.approx(forces, rel=1E-6) == f


def test_nonbonded_kernel():
    pairs = list(itertools.combinations(list(range(top.n_atom)), 2))
    parameters = [(i, j / 1000, i / 10 + j / 10) for i, j in pairs]

    kernel = NonbondedKernel(top.positions,
                             pairs,
                             parameters)
    r, energy, forces = kernel.evaluate()
    print()
    print(r)
    print(energy)
    print(sum(energy))
    print(forces)

    force = mm.NonbondedForce()
    force.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
    force.setForceGroup(5)
    for i in range(top.n_atom):
        force.addParticle(0, 1, 0)
    for (i, j), parameter in zip(pairs, parameters):
        force.addException(i, j, parameter[2], parameter[1], parameter[0])
    e, f = calc_energy_with_omm(force, top.positions)

    assert pytest.approx(sum(energy), rel=1E-6) == e
    assert pytest.approx(forces, rel=1E-6) == f
