#!/usr/bin/env python3

import sys
import random
import openmm.openmm as mm
from openmm import app
import mstools.ommhelper as oh
from mstools.ommhelper.unit import *


def run_simulation(nstep, gro_file='conf.gro', psf_file='topol.psf', prm_file='ff.prm',
                   dt=0.001, T=333, thole=0, restart=None):
    print('Building system...')
    gro = oh.GroFile(gro_file)
    lz = gro.getUnitCellDimensions()[2].value_in_unit(nm)
    psf = oh.OplsPsfFile(psf_file, periodicBoxVectors=gro.getPeriodicBoxVectors())
    prm = app.CharmmParameterSet(prm_file)
    system = psf.createSystem(prm, nonbondedMethod=app.PME, nonbondedCutoff=1.2 * nm,
                              constraints=app.HBonds, rigidWater=True, verbose=True)
    nbforce = next(f for f in system.getForces() if type(f) == mm.NonbondedForce)
    is_drude = any(type(f) == mm.DrudeForce for f in system.getForces())

    ### assign groups
    atoms = list(psf.topology.atoms())
    group_mos = [atom.index for atom in atoms if atom.residue.name == 'MoS2']
    group_mos_core = [i for i in group_mos if not atoms[i].name.startswith('D')]
    group_img = [atom.index for atom in atoms if atom.residue.name == 'IMG']
    group_ils = [atom.index for atom in atoms if atom.residue.name not in ['MoS2', 'IMG']]
    group_ils_not_H = [i for i in group_ils if not atoms[i].name.startswith('H')]
    print('Assigning groups...')
    print('    Number of atoms in group %10s: %i' % ('mos', len(group_mos)))
    print('    Number of atoms in group %10s: %i' % ('ils', len(group_ils)))
    print('    Number of atoms in group %10s: %i' % ('img', len(group_img)))
    print('    Number of atoms in group %10s: %i' % ('mos_core', len(group_mos_core)))
    print('    Number of atoms in group %10s: %i' % ('ils_not_H', len(group_ils_not_H)))

    ### apply TT damping for CLPol force field
    donors = [atom.idx for atom in psf.atom_list if atom.attype == 'HO']
    if is_drude and len(donors) > 0:
        print('Add TT damping between HO and Drude dipoles')
        ttforce = oh.CLPolCoulTT(system, donors)
        print(ttforce.getEnergyFunction())

    ### apply Thole screening between ILs and MoS2
    if thole != 0:
        print('Add Thole screening between ILs and MoS2...')
        parent_drude_dict = dict(psf.drudepair_list)
        drude_parent_dict = {pair[1]: pair[0] for pair in psf.drudepair_list}
        tforce = mm.CustomNonbondedForce('-%.6f*q1*q2/r*(1+0.5*screen*r)*exp(-screen*r);'
                                         'screen=%.4f*(alpha1*alpha2)^(-1/6)' % (
                                             oh.CONST.ONE_4PI_EPS0, thole))
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
        tforce.addInteractionGroup(set(group_mos), set(group_ils_not_H))
        tforce.setForceGroup(9)
        system.addForce(tforce)

    ### apply slab correction, need to move all atoms to primary cell
    print('Add dipole slab correction...')
    for i in range(system.getNumParticles()):
        gro.positions[i] += (0, 0, 1) * nm
    oh.slab_correction(system)

    ### restrain the movement of MoS2 cores (Drude particles are free if exist)
    print('Add restraint for MoS2 cores...')
    oh.spring_self(system, gro.positions.value_in_unit(nm), group_mos_core,
                   [0.01, 0.01, 5.0] * kcal_mol / angstrom ** 2)

    ### velocity-Verlet-middle integrator
    # from velocityverletplugin import VVIntegrator
    # integrator = VVIntegrator(T * kelvin, 10 / ps, 1 * kelvin, 40 / ps, dt * ps)
    # integrator.setUseMiddleScheme(True)
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
    sim.reporters.append(oh.StateDataReporter(sys.stdout, 10000, append=append, box=False))
    sim.reporters.append(oh.DrudeTemperatureReporter('T_drude.txt', 100000, append=append))

    print('Running...')
    oh.energy_decomposition(sim)
    sim.runForClockTime(31.9, 'rst.cpt', 'rst.xml', 1)
    print('# clock time limit: step= %i time= %f' % (
        sim.currentStep, sim.context.getState().getTime().value_in_unit(ps)))


if __name__ == '__main__':
    oh.print_omm_info()
    run_simulation(nstep=int(1E8), T=333, thole=0, restart=None)
