import sys
import simtk.openmm as mm
from simtk.openmm import app
from simtk.unit import kelvin, bar
from simtk.unit import picosecond as ps, nanometer as nm, kilojoule_per_mole as kJ_mol
from mstools.omm import OplsPsfFile, GroFile
from mstools.omm import GroReporter, XmlStateReporter, DrudeTemperatureReporter
from mstools.omm.utils import print_omm_info, minimize, apply_mc_barostat

def run_simulation(nstep, gro_file='conf.gro', psf_file='topol.psf', prm_file='ff.prm',
                   T=300, P=1, tcoupl='langevin', pcoupl='iso'):
    print('Building system...')
    gro = app.GromacsGroFile(gro_file)
    psf = OplsPsfFile(psf_file, periodicBoxVectors=gro.getPeriodicBoxVectors())
    prm = app.CharmmParameterSet(prm_file)
    system = psf.createSystem(prm, nonbondedMethod=app.PME, ewaldErrorTolerance=1E-5, nonbondedCutoff=1.2 * nm,
                              constraints=app.HBonds, rigidWater=True, verbose=True)
    is_drude = any([isinstance(f, mm.DrudeForce) for f in system.getForces()])

    print('Initializing simulation...')
    if tcoupl == 'langevin':
        if is_drude:
            print('    Drude Langevin thermostat: 5.0 /ps, 20 /ps')
            integrator = mm.DrudeLangevinIntegrator(T * kelvin, 5.0 / ps, 1 * kelvin, 20 / ps, 0.001 * ps)
            integrator.setMaxDrudeDistance(0.02 * nm)
        else:
            print('    Langevin thermostat: 1.0 /ps')
            integrator = mm.LangevinIntegrator(T * kelvin, 1.0 / ps, 0.001 * ps)
    elif tcoupl == 'nose-hoover':
        if is_drude:
            print('    Drude Temperature-Grouped Nose-Hoover thermostat: 0.1 ps, 0.025 ps')
            from velocityverletplugin import VVIntegrator
            integrator = VVIntegrator(T * kelvin, 0.1 * ps, 1 * kelvin, 0.025 * ps, 0.001 * ps, 1, 3, True, True)
            integrator.setMaxDrudeDistance(0.02 * nm)
        else:
            print('    Nose-Hoover thermostat: 10 /ps')
            integrator = mm.NoseHooverIntegrator(T * kelvin, 10 / ps, 0.001 * ps)
    else:
        raise Exception('Available thermostat: langevin, nose-hoover')

    if pcoupl is not None:
        apply_mc_barostat(system, pcoupl, P, T)

    _platform = mm.Platform.getPlatformByName('CUDA')
    _properties = {'CudaPrecision': 'mixed'}
    sim = app.Simulation(psf.topology, system, integrator, _platform, _properties)
    sim.context.setPositions(gro.positions)
    sim.context.setVelocitiesToTemperature(T * kelvin)
    sim.reporters.append(XmlStateReporter('state.xml', max(nstep // 10, 100000)))
    sim.reporters.append(GroReporter('dump.gro', 'logfreq', enforcePeriodicBox=False))
    sim.reporters.append(app.DCDReporter('dump.dcd', 10000, enforcePeriodicBox=False))
    sim.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True,
                                               volume=True, density=True, speed=True, separator='\t'))
    if is_drude:
        sim.reporters.append(DrudeTemperatureReporter('T_drude.txt', 10000))

    state = sim.context.getState(getEnergy=True)
    print('Initial Energy: ' + str(state.getPotentialEnergy()))
    print('Minimizing...')
    minimize(sim, 500, gro_out='em.gro')

    print('Running...')
    sim.step(nstep)
    sim.saveCheckpoint('rst.cpt')


if __name__ == '__main__':
    print_omm_info()
    run_simulation(nstep=100000, T=300, P=1, tcoupl='langevin', pcoupl='iso')
