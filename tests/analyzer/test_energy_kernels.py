#!/usr/bin/env python3

import os
import itertools
import tempfile
import filecmp
import pytest
from mstools.topology import Topology, Molecule
from mstools.analyzer.energy_kernels import *

from simtk import openmm as mm

cwd = os.path.dirname(os.path.abspath(__file__))

top = Topology.open(cwd + '/../topology/files/CH3NH2.pdb')
system = mm.System()
for i in range(top.n_atom):
    system.addParticle(1)

force = mm.HarmonicBondForce()
force.setForceGroup(1)
for i, bond in enumerate(top.bonds):
    force.addBond(bond.atom1.id, bond.atom2.id, 0.1 + 0.01 * i, 100 + 10 * i)
system.addForce(force)

force = mm.HarmonicAngleForce()
force.setForceGroup(2)
for i, angle in enumerate(top.angles):
    force.addAngle(angle.atom1.id, angle.atom2.id, angle.atom3.id, (120 + 12 * i) / 180 * np.pi, 100 + 10 * i)
system.addForce(force)

force = mm.PeriodicTorsionForce()
force.setForceGroup(3)
for i, dihedral in enumerate(top.dihedrals):
    for k in range(1, 5):
        force.addTorsion(dihedral.atom1.id, dihedral.atom2.id, dihedral.atom3.id, dihedral.atom4.id,
                         k, ((k - 1) % 2) * np.pi, k + 0.1 * i)
system.addForce(force)

force = mm.CustomTorsionForce('0.5*k*min(dtheta, 2*pi-dtheta)^2;'
                              'dtheta = abs(theta-theta0);'
                              'pi = 3.1415926535')
force.addPerTorsionParameter('theta0')
force.addPerTorsionParameter('k')
force.setForceGroup(4)
for i, dihedral in enumerate(top.dihedrals):
    force.addTorsion(dihedral.atom1.id, dihedral.atom2.id, dihedral.atom3.id, dihedral.atom4.id,
                     # [1 + 0.1 * i, 2 + 0.1 * i]
                     [top.dihedrals[i].evaluate() * 0.8 / 180 * np.pi, 2 + 0.1 * i]
                     )
system.addForce(force)

force = mm.NonbondedForce()
force.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
force.setForceGroup(5)
for i in range(top.n_atom):
    force.addParticle(0, 1, 0)
pairs = list(itertools.combinations(list(range(top.n_atom)), 2))
parameters = [(i, j / 1000, i / 10 + j / 10) for i, j in pairs]
for (i, j), parameter in zip(pairs, parameters):
    force.addException(i, j, parameter[2], parameter[1], parameter[0])
system.addForce(force)

integrator = mm.VerletIntegrator(0.001)
platform = mm.Platform.getPlatformByName('Reference')
context = mm.Context(system, integrator, platform)
context.setPositions(top.positions)


def test_harmonic_bond_kernel():
    kernel = HarmonicBondKernel(top.positions,
                                [(bond.atom1.id, bond.atom2.id) for bond in top.bonds],
                                [[0.1 + 0.01 * i, 100 + 10 * i] for i in range(top.n_bond)],
                                )
    r, energy, forces = kernel.evaluate()
    print()
    print(r)
    print(energy)
    print(sum(energy))
    print(forces)

    state = context.getState(getEnergy=True, getForces=True, groups={1})
    print(state.getPotentialEnergy())
    print(state.getForces(asNumpy=True))


def test_harmonic_angle_kernel():
    kernel = HarmonicAngleKernel(top.positions,
                                 [(angle.atom1.id, angle.atom2.id, angle.atom3.id) for angle in top.angles],
                                 [[120 + 12 * i, 100 + 10 * i] for i in range(top.n_angle)],
                                 )
    theta, energy, forces = kernel.evaluate()
    print()
    print(theta * 180 / np.pi)
    print(energy)
    print(sum(energy))
    print(forces)

    state = context.getState(getEnergy=True, getForces=True, groups={2})
    print(state.getPotentialEnergy())
    print(state.getForces(asNumpy=True))


def test_periodic_torsion_kernel():
    pos = top.positions + np.random.random(top.positions.shape)

    kernel = PeriodicTorsionKernel(pos,
                                   [(dihedral.atom1.id, dihedral.atom2.id, dihedral.atom3.id, dihedral.atom4.id)
                                    for dihedral in top.dihedrals],
                                   [[1 + 0.1 * i, 2 + 0.1 * i, 3 + 0.1 * i, 4 + 0.1 * i]
                                    # [[1 + 0.1 * i, 0, 0, 0]
                                    for i in range(top.n_dihedral)],
                                   )
    phi, energy, forces = kernel.evaluate()
    print()
    for d, v in zip(top.dihedrals, phi * 180 / np.pi):
        print(d.atom1.id + 1, d.atom2.id + 1, d.atom3.id + 1, d.atom4.id + 1, v)
    print(energy)
    print(sum(energy))
    print(forces)

    context.setPositions(pos)
    state = context.getState(getEnergy=True, getForces=True, groups={3})
    print(state.getPotentialEnergy())
    print(state.getForces(asNumpy=True))


def test_harmonic_torsion_kernel():
    pos = top.positions + np.random.random(top.positions.shape)
    print([[top.dihedrals[i].evaluate(), 2 + 0.1 * i]
                                    for i in range(top.n_dihedral)])

    kernel = HarmonicTorsionKernel(pos,
                                   [(dihedral.atom1.id, dihedral.atom2.id, dihedral.atom3.id, dihedral.atom4.id)
                                    for dihedral in top.dihedrals],
                                   # [[1 + 0.1 * i, 2 + 0.1 * i] for i in range(top.n_dihedral)],
                                   [[top.dihedrals[i].evaluate()*0.8, 2 + 0.1 * i]
                                    for i in range(top.n_dihedral)]
                                   )
    phi, energy, forces = kernel.evaluate()
    print()
    print(phi)
    print(energy)
    print(sum(energy))
    print(forces)

    context.setPositions(pos)
    state = context.getState(getEnergy=True, getForces=True, groups={4})
    print(state.getPotentialEnergy())
    print(state.getForces(asNumpy=True))


def test_nonbonded_kernel():
    pairs = list(itertools.combinations(list(range(top.n_atom)), 2))
    parameters = [(i, j / 1000, i / 10 + j / 10) for i, j in pairs]

    kernel = NonbondedKernel(top.positions, pairs, parameters)
    r, energy, forces = kernel.evaluate()
    print()
    print(r)
    print(energy)
    print(sum(energy))
    print(forces)

    state = context.getState(getEnergy=True, getForces=True, groups={5})
    print(state.getPotentialEnergy())
    print(state.getForces(asNumpy=True))
