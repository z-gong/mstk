import sys
import simtk.openmm as mm
from simtk.openmm import app
from simtk.unit import kelvin, bar
from simtk.unit import picosecond as ps, nanometer as nm, kilojoule_per_mole as kJ_mol
from ..omm import OplsPsfFile, GroReporter, GroFile


def gen_simulation(gro_file, psf_file, prm_file, T, P, pcoupl='iso'):
    print('Building system...')
    gro = app.GromacsGroFile(gro_file)
    psf = OplsPsfFile(psf_file, periodicBoxVectors=gro.getPeriodicBoxVectors())
    prm = app.CharmmParameterSet(prm_file)
    system = psf.createSystem(prm, nonbondedMethod=app.PME, ewaldErrorTolerance=1E-5,
                              nonbondedCutoff=1.2 * nm,
                              constraints=app.HBonds, rigidWater=True,
                              verbose=True)

    if not psf.is_drude:
        print('\tLangevin thermostat: 1.0/ps')
        integrator = mm.LangevinIntegrator(T * kelvin, 1.0 / ps, 0.001 * ps)
    else:
        print('\tDual Langevin thermostat: 10/ps, 50/ps')
        integrator = mm.DrudeLangevinIntegrator(T * kelvin, 10 / ps, 1 * kelvin, 50 / ps, 0.001 * ps)
        integrator.setMaxDrudeDistance(0.02 * nm)

    if pcoupl == 'iso':
        print('\tIsotropic barostat')
        system.addForce(mm.MonteCarloBarostat(P * bar, T * kelvin, 25))
    elif pcoupl == 'xyz':
        print('\tAnisotropic barostat')
        system.addForce(mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, True, True, True, 25))
    elif pcoupl == 'xy':
        print('\tAnisotropic Barostat only for X and Y')
        system.addForce(mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, True, True, False, 25))
    else:
        print('\tNo barostat')

    print('Initializing simulation...')
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
    sim = app.Simulation(psf.topology, system, integrator, platform, properties)
    sim.context.setPositions(gro.positions)
    sim.context.setVelocitiesToTemperature(T * kelvin)
    sim.reporters.append(GroReporter('dump.gro', 10000))
    sim.reporters.append(app.DCDReporter('dump.dcd', 1000))
    sim.reporters.append(app.StateDataReporter(sys.stdout, 500, step=True, temperature=True,
                                               potentialEnergy=True, kineticEnergy=True,
                                               volume=True, density=True,
                                               elapsedTime=True, speed=True, separator='\t'))
    state = sim.context.getState(getEnergy=True)
    print('Energy: ' + str(state.getPotentialEnergy()))
    return sim


def print_omm_info():
    print(mm.__version__)
    print(mm.version.openmm_library_path)
    print(mm.pluginLoadedLibNames)
    print(mm.Platform.getPluginLoadFailures())
    print([mm.Platform.getPlatform(index).getName() for index in range(mm.Platform.getNumPlatforms())])


def minimize(sim, tolerance):
    print('Minimizing...')
    sim.minimizeEnergy(tolerance=tolerance * kJ_mol)
    state = sim.context.getState(getPositions=True, getEnergy=True)
    with open('em.gro', 'w') as f:
        GroFile.writeFile(sim.topology, state.getTime(),
                          state.getPositions(), state.getPeriodicBoxVectors(), f)
    print('Energy: ' + str(state.getPotentialEnergy()))


def run_simulation(nstep, gro='conf.gro', psf='topol.psf', prm='ff.prm',
                   T=300, P=1, pcoupl='iso'):
    print_omm_info()
    sim = gen_simulation(gro, psf, prm, T, P, pcoupl)
    minimize(sim, 200)
    print('Running...')
    sim.step(nstep)


if __name__ == '__main__':
    run_simulation(nstep=10000)
