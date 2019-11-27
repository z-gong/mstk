import sys
import simtk.openmm as mm
from simtk.openmm import app
from simtk.unit import kelvin, bar, elementary_charge as e
from simtk.unit import picosecond as ps, nanometer as nm, kilojoule_per_mole as kJ_mol
from mstools.omm import OplsPsfFile, GroReporter, GroFile
from mstools.omm.drudetemperaturereporter import DrudeTemperatureReporter


def slab_correction(system: mm.System):
    muz = mm.CustomExternalForce('q*z')
    muz.addPerParticleParameter('q')
    nbforce = [system.getForce(i) for i in range(system.getNumForces())
               if system.getForce(i).__class__ == mm.NonbondedForce][0]
    qsum = 0
    for i in range(nbforce.getNumParticles()):
        q = nbforce.getParticleParameters(i)[0].value_in_unit(e)
        muz.addParticle(i, [q])
        qsum += q
    if abs(qsum) > 1E-4:
        raise Exception('Slab correction is not valid for charged system')

    box = system.getDefaultPeriodicBoxVectors()
    vol = (box[0][0] * box[1][1] * box[2][2]).value_in_unit(nm ** 3)
    cvforce = mm.CustomCVForce('k*muz*muz')
    # convert from e^2/nm to kJ/mol
    _conv = (1.602E-19) ** 2 / 1E-9 / (4 * 3.1415926) / 8.854E-12 * 6.022E23 / 1000
    cvforce.addGlobalParameter('k', 2 * 3.1415926 / vol * _conv)
    cvforce.addCollectiveVariable('muz', muz)
    system.addForce(cvforce)


def gen_simulation(gro_file, psf_file, prm_file, T, P, tcoupl='langevin', pcoupl='iso', slab=False, platform='CUDA'):
    print('Building system...')
    gro = app.GromacsGroFile(gro_file)
    psf = OplsPsfFile(psf_file, periodicBoxVectors=gro.getPeriodicBoxVectors())
    prm = app.CharmmParameterSet(prm_file)
    system = psf.createSystem(prm, nonbondedMethod=app.PME, ewaldErrorTolerance=1E-5,
                              nonbondedCutoff=1.2 * nm,
                              constraints=app.HBonds, rigidWater=True,
                              verbose=True)

    print('Initializing simulation...')
    if tcoupl == 'langevin':
        if psf.is_drude:
            print('\tDrude Langevin thermostat: 5.0/ps, 20/ps')
            integrator = mm.DrudeLangevinIntegrator(T * kelvin, 5.0 / ps, 1 * kelvin, 20 / ps, 0.001 * ps)
        else:
            print('\tLangevin thermostat: 1.0/ps')
            integrator = mm.LangevinIntegrator(T * kelvin, 1.0 / ps, 0.001 * ps)
    elif tcoupl == 'nose-hoover':
        if psf.is_drude:
            print('\tDrude Nose-Hoover thermostat: 10/ps, 40/ps')
            integrator = mm.DrudeNoseHooverIntegrator(T * kelvin, 10 / ps, 1 * kelvin, 40 / ps, 0.001 * ps)
        else:
            print('\tNose-Hoover thermostat: 10/ps')
            integrator = mm.NoseHooverIntegrator(T * kelvin, 10 / ps, 0.001 * ps)
    else:
        raise Exception('Available thermostat: langevin, nose-hoover')
    integrator.setMaxDrudeDistance(0.025 * nm)

    if pcoupl == 'iso':
        print('\tIsotropic barostat')
        system.addForce(mm.MonteCarloBarostat(P * bar, T * kelvin, 25))
    elif pcoupl == 'xyz':
        print('\tAnisotropic barostat')
        system.addForce(mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, True, True, True, 25))
    elif pcoupl == 'xy':
        print('\tAnisotropic barostat only for X and Y')
        system.addForce(mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, True, True, False, 25))
    elif pcoupl == 'z':
        print('\tAnisotropic barostat only for Z')
        system.addForce(mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, False, False, True, 25))
    elif pcoupl is None:
        print('\tNo barostat')
    else:
        raise Exception('Available barostat: iso, xyz, xy, z, None')

    if slab:
        print('\tSlab correction applied')
        slab_correction(system)

    properties = {}
    if platform == 'CUDA':
        print('\tPlatform: CUDA')
        _platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed'}
    elif platform == 'OpenCL':
        print('\tPlatform: OpenCL')
        _platform = mm.Platform.getPlatformByName('OpenCL')
        properties = {'OpenCLPrecision': 'mixed'}
    elif platform == 'Reference':
        print('\tPlatform: Reference')
        _platform = mm.Platform.getPlatformByName('Reference')
    else:
        raise Exception('Available platform: CUDA, OpenCL, Reference')

    sim = app.Simulation(psf.topology, system, integrator, _platform, properties)
    sim.context.setPositions(gro.positions)
    sim.context.setVelocitiesToTemperature(T * kelvin)
    sim.reporters.append(GroReporter('dump.gro', 10000))
    sim.reporters.append(app.DCDReporter('dump.dcd', 1000, enforcePeriodicBox=False))
    sim.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, temperature=True,
                                               potentialEnergy=True, kineticEnergy=True,
                                               volume=True, density=True,
                                               elapsedTime=True, speed=True, separator='\t'))
    if psf.is_drude:
        sim.reporters.append(DrudeTemperatureReporter('T_drude.txt', 10000))

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
                   T=300, P=1, tcoupl='langevin', pcoupl='iso', slab=False,
                   platform='OpenCL'):
    print_omm_info()
    sim = gen_simulation(gro, psf, prm, T, P, tcoupl, pcoupl, slab, platform)
    minimize(sim, 200)
    print('Running...')
    sim.step(nstep)
    sim.saveCheckpoint('rst.cpt')
    sim.saveState('rst.xml')


if __name__ == '__main__':
    run_simulation(nstep=10000)
