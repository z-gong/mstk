#!/usr/bin/env python3

import sys
import random
import openmm.openmm as mm
from openmm.openmm import app
import mstk.ommhelper as oh
from mstk.ommhelper.unit import *


def run_simulation(nstep, gro_file='conf.gro', psf_file='topol.psf', prm_file='ff.prm',
                   dt=0.001, T=333, voltage=0, screen=0, restart=None):
    print('Building system...')
    gro = oh.GroFile(gro_file)
    lz = gro.getUnitCellDimensions()[2].value_in_unit(nm)
    psf = oh.OplsPsfFile(psf_file, periodicBoxVectors=gro.getPeriodicBoxVectors())
    prm = app.CharmmParameterSet(prm_file)
    system = psf.createSystem(prm, nonbondedMethod=app.PME, nonbondedCutoff=1.2 * nm,
                              constraints=app.HBonds, rigidWater=True, verbose=True)
    nbforce = next(f for f in system.getForces() if type(f) == mm.NonbondedForce)
    ljforce = next(f for f in system.getForces() if type(f) == mm.CustomNonbondedForce)
    bforce = next(f for f in system.getForces() if type(f) == mm.HarmonicBondForce)
    is_drude = any(type(f) == mm.DrudeForce for f in system.getForces())

    ### assign groups
    atoms = list(psf.topology.atoms())
    group_mos = [atom.index for atom in atoms if atom.residue.name == 'MoS2']
    group_mos_core = [i for i in group_mos if not atoms[i].name.startswith('D')]
    group_img = [atom.index for atom in atoms if atom.residue.name == 'IMG']
    group_ils = [atom.index for atom in atoms if atom.residue.name not in ['MoS2', 'IMG']]
    group_ils_drude = [i for i in group_ils if atoms[i].name.startswith('D')]
    image_pairs = list(zip(group_ils, group_img))
    print('Assigning groups...')
    print('    Number of atoms in group %10s: %i' % ('mos', len(group_mos)))
    print('    Number of atoms in group %10s: %i' % ('ils', len(group_ils)))
    print('    Number of atoms in group %10s: %i' % ('img', len(group_img)))
    print('    Number of atoms in group %10s: %i' % ('mos_core', len(group_mos_core)))

    ### apply TT damping for CLPol force field
    donors = [atom.idx for atom in psf.atom_list if atom.attype == 'HO']
    if is_drude and len(donors) > 0:
        print('Add TT damping between HO and Drude dipoles')
        ttforce = oh.CLPolCoulTT(system, donors)
        print(ttforce.getEnergyFunction())

    ### assign image charges
    for parent, image in image_pairs:
        q, _, _ = nbforce.getParticleParameters(parent)
        nbforce.setParticleParameters(image, -q, 1, 0)

    ### eliminate the interactions between images and MoS2
    ljforce.addInteractionGroup(set(group_img), set(group_ils))
    ljforce.addInteractionGroup(set(group_mos + group_ils), set(group_mos + group_ils))

    ### apply screening between ils and images
    if screen != 0:
        print('Add screening between ILs and images...')
        tforce = mm.CustomNonbondedForce('-%.6f*q1*q2/r*(1+0.5*screen*r)*exp(-screen*r);'
                                         'screen=%.4f' % (oh.CONST.ONE_4PI_EPS0, screen))
        print(tforce.getEnergyFunction())
        tforce.addPerParticleParameter('q')
        for i in range(system.getNumParticles()):
            q, _, _ = nbforce.getParticleParameters(i)
            tforce.addParticle([q])
        tforce.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        tforce.setCutoffDistance(1.2 * nm)
        tforce.addInteractionGroup(set(group_ils), set(group_img))
        tforce.setForceGroup(9)
        system.addForce(tforce)

    ### restrain the movement of MoS2 cores (Drude particles are free if exist)
    print('Add restraint for MoS2...')
    oh.spring_self(system, gro.positions.value_in_unit(nm), group_mos_core,
                   [0.01, 0.01, 5.0] * kcal_mol / angstrom ** 2)

    ### in case Drude particles kiss their images
    print('Add wall for Drude particles of ILs...')
    wall = oh.wall_lj126(system, group_ils_drude, 'z', [0, lz / 2],
                         epsilon=0.5 * kcal_mol, sigma=0.15 * nm)
    print(wall.getEnergyFunction())

    ### randomlize particle positions to get ride of overlap
    random.seed(0)
    for i in range(system.getNumParticles()):
        gro.positions[i] += (random.random(), random.random(), random.random()) * nm / 1000

    ### velocity-Verlet-middle integrator
    from velocityverletplugin import VVIntegrator
    integrator = VVIntegrator(T * kelvin, 10 / ps, 1 * kelvin, 40 / ps, dt * ps)
    integrator.setUseMiddleScheme(True)
    integrator.setMaxDrudeDistance(0.02 * nm)
    ### thermostat MoS2 by Langevin dynamics
    for i in group_mos:
        integrator.addParticleLangevin(i)
    ### assign image pairs
    integrator.setMirrorLocation(lz / 2 * nm)
    for parent, image in image_pairs:
        integrator.addImagePair(image, parent)
        # add fake bond between image and parent so that they are always in the same periodic cell
        bforce.addBond(image, parent, 0, 0)
    ### apply electric field on ils
    if voltage != 0:
        integrator.setElectricField(voltage / lz * 2 * volt / nm)
        for i in group_ils:
            integrator.addParticleElectrolyte(i)

    print('Initializing simulation...')
    _platform = mm.Platform.getPlatformByName('CUDA')
    _properties = {'CudaPrecision': 'mixed'}
    sim = app.Simulation(psf.topology, system, integrator, _platform, _properties)
    if restart:
        sim.loadCheckpoint(restart)
        sim.currentStep = round(sim.context.getState().getTime().value_in_unit(ps) / dt / 10) * 10
        sim.context.setTime(sim.currentStep * dt)
        append = True
    else:
        sim.context.setPositions(gro.positions)
        sim.context.setVelocitiesToTemperature(T * kelvin)
        append = False

    sim.reporters.append(app.DCDReporter('dump.dcd', 10000, enforcePeriodicBox=False,
                                         append=append))
    sim.reporters.append(oh.CheckpointReporter('cpt.cpt', 10000))
    sim.reporters.append(oh.GroReporter('dump.gro', 1000, logarithm=True,
                                        subset=group_mos + group_ils, append=append))
    sim.reporters.append(oh.StateDataReporter(sys.stdout, 10000, box=False, append=append))
    sim.reporters.append(oh.DrudeTemperatureReporter('T_drude.txt', 100000, append=append))

    print('Running...')
    oh.energy_decomposition(sim)
    sim.runForClockTime(99.9, 'rst.cpt', 'rst.xml', 1)
    print('# clock time limit: step= %i time= %f' % (
        sim.currentStep, sim.context.getState().getTime().value_in_unit(ps)))


if __name__ == '__main__':
    oh.print_omm_info()
    run_simulation(nstep=int(1E8), T=333, voltage=5, screen=0, restart=None)
