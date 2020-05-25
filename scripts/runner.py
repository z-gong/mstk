#!/usr/bin/env python3

import sys
import signal
import simtk.openmm as mm
from simtk.openmm import app
import mstools.ommhelper as oh
from mstools.ommhelper.unit import *


class SignalHandler():
    def __init__(self):
        self.SIGINT = False

    def sigint_handler(self, signal, frame):
        self.SIGINT = True


sig_handler = SignalHandler()
signal.signal(signal.SIGINT, sig_handler.sigint_handler)
signal.signal(signal.SIGTERM, sig_handler.sigint_handler)


def run_simulation(nstep, gro_file='conf.gro', psf_file='topol.psf', prm_file='ff.prm',
                   dt=0.001, T=300, P=1, tcoupl='langevin', pcoupl='iso', restart=None):
    print('Building system...')
    gro = oh.GroFile(gro_file)
    psf = oh.OplsPsfFile(psf_file, periodicBoxVectors=gro.getPeriodicBoxVectors())
    prm = app.CharmmParameterSet(prm_file)
    system = psf.createSystem(prm, nonbondedMethod=app.PME, nonbondedCutoff=1.2 * nm,
                              constraints=app.HBonds, rigidWater=True, verbose=True)
    is_drude = any(type(f) == mm.DrudeForce for f in system.getForces())

    print('Initializing simulation...')
    if tcoupl == 'langevin':
        if is_drude:
            print('Drude Langevin thermostat: 5.0 /ps, 20 /ps')
            integrator = mm.DrudeLangevinIntegrator(T * kelvin, 5.0 / ps, 1 * kelvin, 20 / ps, dt * ps)
            integrator.setMaxDrudeDistance(0.02 * nm)
        else:
            print('Langevin thermostat: 1.0 /ps')
            integrator = mm.LangevinIntegrator(T * kelvin, 1.0 / ps, dt * ps)
    elif tcoupl == 'nose-hoover':
        if is_drude:
            print('Drude temperature-grouped Nose-Hoover thermostat: 10 /ps, 40 /ps')
        else:
            print('Nose-Hoover thermostat: 10 /ps')
        from velocityverletplugin import VVIntegrator
        integrator = VVIntegrator(T * kelvin, 10 / ps, 1 * kelvin, 40 / ps, dt * ps)
        integrator.setMaxDrudeDistance(0.02 * nm)
    else:
        raise Exception('Available thermostat: langevin, nose-hoover')

    if pcoupl is not None:
        oh.apply_mc_barostat(system, pcoupl, P, T)

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
        oh.energy_decomposition(sim)
        # minimize(sim, 100, 'em.gro')
        append = False

    sim.reporters.append(app.DCDReporter('dump.dcd', 10000, enforcePeriodicBox=False, append=append))
    sim.reporters.append(oh.CheckpointReporter('cpt.cpt', 10000))
    sim.reporters.append(oh.GroReporter('dump.gro', 'logfreq', append=append))
    sim.reporters.append(oh.StateDataReporter(sys.stdout, 1000, box=False, volume=True, append=append))
    if is_drude:
        sim.reporters.append(oh.DrudeTemperatureReporter('T_drude.txt', 10000, append=append))

    print('Running...')
    i_step = 0
    while i_step < nstep:
        if sig_handler.SIGINT:
            sim.saveCheckpoint('rst.cpt')
            print('# SIGINT step= %i time= %f' % (
                sim.currentStep, sim.context.getState().getTime().value_in_unit(ps)))
            break
        sim.step(10)
        i_step += 10
    sim.saveCheckpoint('rst.cpt')
    sim.saveState('rst.xml')


if __name__ == '__main__':
    oh.print_omm_info()
    run_simulation(nstep=100000, T=300, P=1, tcoupl='langevin', pcoupl='iso', restart=None)
