#!/usr/bin/env python3

import sys
import random
import simtk.openmm as mm
from simtk.openmm import app
from simtk.unit import kelvin, bar, volt
from simtk.unit import picosecond as ps, nanometer as nm, kilojoule_per_mole as kJ_mol
from mstools.omm import OplsPsfFile, GroReporter, GroFile
from mstools.omm.drudetemperaturereporter import DrudeTemperatureReporter
from mstools.omm.forces import spring_self, wall_power, wall_lj126
from mstools.omm.utils import print_omm_info, minimize


def run_simulation(nstep, gro_file='conf.gro', psf_file='topol.psf', prm_file='ff.prm', T=333, voltage=0):
    print('Building system...')
    gro = app.GromacsGroFile(gro_file)
    box_z = gro.getUnitCellDimensions()[2].value_in_unit(nm)
    psf = OplsPsfFile(psf_file, periodicBoxVectors=gro.getPeriodicBoxVectors())
    prm = app.CharmmParameterSet(prm_file)
    system = psf.createSystem(prm, nonbondedMethod=app.PME, ewaldErrorTolerance=1E-5,
                              nonbondedCutoff=1.2 * nm, constraints=app.HBonds, rigidWater=True,
                              verbose=True)
    forces = {system.getForce(i).__class__.__name__: system.getForce(i) for i in range(system.getNumForces())}

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

    ### move most of the atoms into first periodic box
    for i in range(system.getNumParticles()):
        gro.positions[i] += (0, 0, box_z / 2) * nm

    ### restrain the movement of MoS2 cores (Drude particles are free if exist)
    spring_self(system, gro.positions.value_in_unit(nm), group_mos_core, 0.001 * 418.4, 0.001 * 418.4, 5.000 * 418.4)
    wall = wall_lj126(system, group_ils_drude, 'z', [box_z / 2, box_z], epsilon=0.5 * 4.184, sigma=0.2)
    # wall = wall_power(system, group_ils_drude, 'z', [box_z / 2 + 0.2, box_z - 0.2], k=0.5 * 4.184, cutoff=0.0245)
    print(wall.getEnergyFunction())

    ### randomlize particle positions to get ride of overlap
    random.seed(0)
    for i in range(system.getNumParticles()):
        gro.positions[i] += (random.random(), random.random(), random.random()) * nm / 1000

    ### eliminate the interactions between images and MoS2
    cforce = forces['CustomNonbondedForce']
    cforce.addInteractionGroup(group_img, group_ils)
    cforce.addInteractionGroup(group_mos + group_ils, group_mos + group_ils)

    ### TGNH thermostat for ils
    from drudenoseplugin import DrudeNoseHooverIntegrator
    integrator = DrudeNoseHooverIntegrator(T * kelvin, 0.1 * ps, 1 * kelvin, 0.025 * ps, 0.001 * ps, 1, 3, True, True)
    integrator.setMaxDrudeDistance(0.02 * nm)
    ### thermostat MoS2 by Langevin dynamics
    for i in group_mos:
        integrator.addParticleLangevin(i)
    ### assign image pairs
    integrator.setMirrorLocation(box_z / 2 * nm)
    for i in group_img:
        integrator.addImagePair(i, i - group_img[0] + group_ils[0])
        # add fake bond between image and parent so that they are always in the same periodic cell
        forces['HarmonicBondForce'].addBond(i, i - group_img[0] + group_ils[0], 1 * nm, 0 * kJ_mol / nm ** 2)
    ### assign image charges
    nbforce = forces['NonbondedForce']
    for i in group_img:
        charge, sig, eps = nbforce.getParticleParameters(i - group_img[0] + group_ils[0])
        nbforce.setParticleParameters(i, -charge, 1 * nm, 0 * kJ_mol)
    ### apply electric field
    integrator.setElectricField(voltage / box_z * 2 * volt / nm)
    for i in group_ils:
        integrator.addParticleElectrolyte(i)
    print('Image charge constant voltage simulation:')
    print('    lz : %s, mirror: %s, electric field: %s' % (
        box_z * nm, integrator.getMirrorLocation(), integrator.getElectricField()))

    print('Initializing simulation...')
    _platform = mm.Platform.getPlatformByName('CUDA')
    _properties = {'CudaPrecision': 'mixed'}
    sim = app.Simulation(psf.topology, system, integrator, _platform, _properties)
    sim.context.setPositions(gro.positions)
    sim.context.setVelocitiesToTemperature(T * kelvin, 12345)
    sim.reporters.append(GroReporter('dump.gro', 100000, enforcePeriodicBox=False))
    sim.reporters.append(app.DCDReporter('dump.dcd', 10000, enforcePeriodicBox=False))
    sim.reporters.append(app.CheckpointReporter('cpt.cpt', 10000))
    sim.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, temperature=True,
                                               potentialEnergy=True, kineticEnergy=True, volume=True, density=True,
                                               elapsedTime=False, speed=True, separator='\t'))
    sim.reporters.append(DrudeTemperatureReporter('T_drude.txt', 10000))

    state = sim.context.getState(getEnergy=True)
    print('Initial Energy: ' + str(state.getPotentialEnergy()))
    # print('Minimizing...')
    # minimize(sim, 500, gro_out='em.gro')
    print('Running...')
    sim.step(nstep)
    sim.saveCheckpoint('rst.cpt')
    sim.saveState('rst.xml')


if __name__ == '__main__':
    print_omm_info()
    run_simulation(nstep=500000, T=333, voltage=5)
