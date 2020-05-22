#!/usr/bin/env python3

import sys
import signal
import random
import simtk.openmm as mm
from simtk.openmm import app
from simtk.unit import kelvin, bar, volt, elementary_charge as qe
from simtk.unit import picosecond as ps, nanometer as nm, angstrom as A
from simtk.unit import kilojoule_per_mole as kJ_mol, kilocalorie_per_mole as kCal_mol
from mstools.omm import OplsPsfFile, GroFile
from mstools.omm import GroReporter, CheckpointReporter, DrudeTemperatureReporter, StateDataReporter
from mstools.omm import spring_self, slab_correction, print_omm_info, minimize, energy_decomposition
from mstools.constant import *


def run_simulation(nstep, gro_file='conf.gro', psf_file='topol.psf', prm_file='ff.prm',
                   dt=0.001, T=333, restart=None):
    print('Building system...')
    gro = GroFile(gro_file)
    lz = gro.getUnitCellDimensions()[2].value_in_unit(nm)
    psf = OplsPsfFile(psf_file, periodicBoxVectors=gro.getPeriodicBoxVectors())
    prm = app.CharmmParameterSet(prm_file)
    system = psf.createSystem(prm, nonbondedMethod=app.PME, nonbondedCutoff=1.2 * nm,
                              constraints=app.HBonds, rigidWater=True, verbose=True)
    forces = {force.__class__.__name__: force for force in system.getForces()}
    nbforce: mm.NonbondedForce = forces['NonbondedForce']
    cforce: mm.CustomNonbondedForce = forces['CustomNonbondedForce']

    ### assign groups
    atoms = list(psf.topology.atoms())
    group_mos = [atom.index for atom in atoms if atom.residue.name == 'MoS2']
    group_mos_core = [atom.index for atom in atoms if
                      atom.residue.name == 'MoS2' and not atom.name.startswith('D')]
    group_img = [atom.index for atom in atoms if atom.residue.name == 'IMG']
    group_ils = [atom.index for atom in atoms if atom.residue.name not in ['MoS2', 'IMG']]
    group_ils_drude = [atom.index for atom in atoms if
                       atom.residue.name not in ['MoS2', 'IMG'] and atom.name.startswith('D')]
    group_ils_not_H = [atom.index for atom in atoms if
                       atom.residue.name not in ['MoS2', 'IMG'] and not atom.name.startswith('H')]
    print('Assigning groups...')
    print('    Number of atoms in group %10s: %i' % ('mos', len(group_mos)))
    print('    Number of atoms in group %10s: %i' % ('ils', len(group_ils)))
    print('    Number of atoms in group %10s: %i' % ('img', len(group_img)))
    print('    Number of atoms in group %10s: %i' % ('mos_core', len(group_mos_core)))
    print('    Number of atoms in group %10s: %i' % ('ils_not_H', len(group_ils_not_H)))

    ### apply thole screening between ILs and MoS2
    print('Add thole screening between ILs and MoS2...')
    thole = 2.6
    parent_drude_dict = dict(psf.drudepair_list)
    drude_parent_dict = {pair[1]: pair[0] for pair in psf.drudepair_list}
    tforce = mm.CustomNonbondedForce('-%.6f*q1*q2/r*(1+0.5*screen*r)*exp(-screen*r);'
                                     'screen=%.4f*(alpha1*alpha2)^(-1/6)' % (ONE_4PI_EPS0, thole))
    print(tforce.getEnergyFunction())
    tforce.addPerParticleParameter('q')
    tforce.addPerParticleParameter('alpha')
    for i in range(system.getNumParticles()):
        if i in parent_drude_dict:
            q, _, _ = nbforce.getParticleParameters(parent_drude_dict[i])
            q = -q
            alpha, _ = psf.drudeconsts_list[i]
        elif i in drude_parent_dict:
            q, _, _ = nbforce.getParticleParameters(i)
            alpha, _ = psf.drudeconsts_list[drude_parent_dict[i]]
        else:
            q = 0 * qe
            alpha = -1
        q = q.value_in_unit(qe)
        alpha = -alpha / 1000  # convert from A^3 to nm^3
        tforce.addParticle([q, alpha])
    tforce.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    tforce.setCutoffDistance(1.2 * nm)
    system.addForce(tforce)
    tforce.addInteractionGroup(set(group_mos), set(group_ils_not_H))

    ### apply slab correction, need to move all atoms to primary cell
    for i in range(system.getNumParticles()):
        gro.positions[i] += (0, 0, 1) * nm
    slab_correction(system)

    ### restrain the movement of MoS2 cores (Drude particles are free if exist)
    spring_self(system, gro.positions.value_in_unit(nm), group_mos_core,
                [0.001, 0.001, 5.0] * kCal_mol / A ** 2)

    ### TGNH thermostat for ils
    # from velocityverletplugin import VVIntegrator
    # integrator = VVIntegrator(T * kelvin, 10 / ps, 1 * kelvin, 40 / ps, dt * ps)
    # integrator.setMaxDrudeDistance(0.02 * nm)
    ### thermostat MoS2 by Langevin dynamics
    # for i in group_mos:
    #     integrator.addParticleLangevin(i)
    ### Langevin thermostat
    integrator = mm.DrudeLangevinIntegrator(T * kelvin, 5 / ps, 1 * kelvin, 20 / ps, dt * ps)
    integrator.setMaxDrudeDistance(0.02 * nm)

    print('Initializing simulation...')
    _platform = mm.Platform.getPlatformByName('CUDA')
    _properties = {'CudaPrecision': 'mixed'}
    sim = app.Simulation(psf.topology, system, integrator, _platform, _properties)
    if restart:
        sim.loadCheckpoint(restart)
        sim.currentStep = int(round(sim.context.getState().getTime().value_in_unit(ps) / dt))
        append = True
    else:
        sim.context.setPositions(gro.positions)
        sim.context.setVelocitiesToTemperature(T * kelvin)
        energy_decomposition(sim)
        # minimize(sim, 100, 'em.gro')
        append = False

    sim.reporters.append(
        app.DCDReporter('dump.dcd', 10000, enforcePeriodicBox=False, append=append))
    sim.reporters.append(CheckpointReporter('cpt.cpt', 10000))
    sim.reporters.append(
        GroReporter('dump.gro', 'logfreq', subset=group_mos + group_ils, append=append))
    sim.reporters.append(StateDataReporter(sys.stdout, 10000, append=append, box=False))
    sim.reporters.append(DrudeTemperatureReporter('T_drude.txt', 100000, append=append))

    print('Running...')
    sim.runForClockTime(31.9, 'rst.cpt', 'rst.xml', 1)


if __name__ == '__main__':
    print_omm_info()
    run_simulation(nstep=int(1E8), T=333, restart=None)
