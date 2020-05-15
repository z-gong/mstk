#!/usr/bin/env python3

import sys
import random
import simtk.openmm as mm
from simtk.openmm import app
from simtk.unit import kelvin, bar, volt
from simtk.unit import picosecond as ps, nanometer as nm, angstrom as A
from simtk.unit import kilojoule_per_mole as kJ_mol, kilocalorie_per_mole as kCal_mol
from mstools.omm import OplsPsfFile, GroFile
from mstools.omm import GroReporter, CheckpointReporter, DrudeTemperatureReporter, StateDataReporter
from mstools.omm import spring_self, wall_lj126, print_omm_info, minimize, energy_decomposition


def run_simulation(nstep, gro_file='conf.gro', psf_file='topol.psf', prm_file='ff.prm',
                   T=333, voltage=0):
    print('Building system...')
    gro = GroFile(gro_file)
    lz = gro.getUnitCellDimensions()[2].value_in_unit(nm)
    psf = OplsPsfFile(psf_file, periodicBoxVectors=gro.getPeriodicBoxVectors())
    prm = app.CharmmParameterSet(prm_file)
    system = psf.createSystem(prm, nonbondedMethod=app.PME, nonbondedCutoff=1.2 * nm,
                              constraints=app.HBonds, rigidWater=True, verbose=True)
    forces = {force.__class__.__name__: force for force in system.getForces()}
    nbforce = forces['NonbondedForce']
    cforce = forces['CustomNonbondedForce']

    ### assign groups
    atoms = list(psf.topology.atoms())
    group_mos = [atom.index for atom in atoms if atom.residue.name == 'MoS2']
    group_mos_core = [atom.index for atom in atoms if
                      atom.residue.name == 'MoS2' and not atom.name.startswith('D')]
    group_img = [atom.index for atom in atoms if atom.residue.name == 'IMG']
    group_ils = [atom.index for atom in atoms if atom.residue.name not in ['MoS2', 'IMG']]
    group_ils_drude = [atom.index for atom in atoms if
                       atom.residue.name not in ['MoS2', 'IMG'] and atom.name.startswith('D')]
    print('Assigning groups...')
    print('    Number of atoms in group %10s: %i' % ('mos', len(group_mos)))
    print('    Number of atoms in group %10s: %i' % ('ils', len(group_ils)))
    print('    Number of atoms in group %10s: %i' % ('img', len(group_img)))
    print('    Number of atoms in group %10s: %i' % ('mos_core', len(group_mos_core)))

    ### restrain the movement of MoS2 cores (Drude particles are free if exist)
    spring_self(system, gro.positions.value_in_unit(nm), group_mos_core,
                [0.001, 0.001, 5.0] * kCal_mol / A ** 2)
    wall = wall_lj126(system, group_ils_drude, 'z', [0, lz / 2],
                      epsilon=0.5 * kCal_mol, sigma=0.15 * nm)
    print(wall.getEnergyFunction())

    ### randomlize particle positions to get ride of overlap
    random.seed(0)
    for i in range(system.getNumParticles()):
        gro.positions[i] += (random.random(), random.random(), random.random()) * nm / 1000

    ### eliminate the interactions between images and MoS2
    cforce.addInteractionGroup(group_img, group_ils)
    cforce.addInteractionGroup(group_mos + group_ils, group_mos + group_ils)

    ### TGNH thermostat for ils
    from velocityverletplugin import VVIntegrator
    integrator = VVIntegrator(T * kelvin, 10 / ps, 1 * kelvin, 40 / ps, 0.001 * ps)
    integrator.setMaxDrudeDistance(0.02 * nm)
    ### thermostat MoS2 by Langevin dynamics
    for i in group_mos:
        integrator.addParticleLangevin(i)
    ### assign image pairs
    integrator.setMirrorLocation(lz / 2 * nm)
    for i in group_img:
        integrator.addImagePair(i, i - group_img[0] + group_ils[0])
        # add fake bond between image and parent so that they are always in the same periodic cell
        forces['HarmonicBondForce'].addBond(i, i - group_img[0] + group_ils[0], 1, 0)
    ### assign image charges
    for i in group_img:
        charge, sig, eps = nbforce.getParticleParameters(i - group_img[0] + group_ils[0])
        nbforce.setParticleParameters(i, -charge, 1, 0)
    ### apply electric field
    if voltage != 0:
        integrator.setElectricField(voltage / lz * 2 * volt / nm)
        for i in group_ils:
            integrator.addParticleElectrolyte(i)

    print('Initializing simulation...')
    _platform = mm.Platform.getPlatformByName('CUDA')
    _properties = {'CudaPrecision': 'mixed'}
    sim = app.Simulation(psf.topology, system, integrator, _platform, _properties)
    sim.context.setPositions(gro.positions)
    sim.context.setVelocitiesToTemperature(T * kelvin)
    sim.reporters.append(app.DCDReporter('dump.dcd', 10000, enforcePeriodicBox=False))
    sim.reporters.append(CheckpointReporter('cpt.cpt', 10000))
    sim.reporters.append(GroReporter('dump.gro', 'logfreq', subset=group_mos + group_ils))
    sim.reporters.append(StateDataReporter(sys.stdout, 10000))
    sim.reporters.append(DrudeTemperatureReporter('T_drude.txt', 100000))

    energy_decomposition(sim)
    minimize(sim, 100, 'em.gro')

    print('Running...')
    sim.step(nstep)
    sim.saveCheckpoint('rst.cpt')
    sim.saveState('rst.xml')
    # sim.runForClockTime(99.9, 'rst.cpt', 'rst.xml', 1)


if __name__ == '__main__':
    print_omm_info()
    run_simulation(nstep=500000, T=333, voltage=5)
