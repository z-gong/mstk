import simtk.openmm as mm
from simtk.unit import kelvin, bar, nanometer as nm, picosecond as ps
from simtk.unit import kilojoule_per_mole as kj_mol, kilocalorie_per_mole as kcal_mol
from .grofile import GroFile


def print_omm_info():
    print(mm.__version__)
    print(mm.version.openmm_library_path)
    print([mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())])
    print(mm.Platform.getPluginLoadFailures())


def minimize(sim, tolerance, gro_out=None):
    sim.minimizeEnergy(tolerance=tolerance * kj_mol)
    state = sim.context.getState(getPositions=True, getEnergy=True)
    print('Minimized energy: ' + str(state.getPotentialEnergy()))

    if gro_out is not None:
        with open(gro_out, 'w') as f:
            GroFile.writeFile(sim.topology, state.getTime(), state.getPositions(),
                              state.getPeriodicBoxVectors(), f)


def apply_mc_barostat(system, pcoupl, P, T, nsteps=100):
    if pcoupl == 'iso':
        print('Isotropic barostat')
        system.addForce(mm.MonteCarloBarostat(P * bar, T * kelvin, nsteps))
    elif pcoupl == 'semi-iso':
        print('Anisotropic barostat with coupled XY')
        system.addForce(mm.MonteCarloMembraneBarostat(P * bar, 0 * bar * nm, T * kelvin,
                                                      mm.MonteCarloMembraneBarostat.XYIsotropic,
                                                      mm.MonteCarloMembraneBarostat.ZFree, nsteps))
    elif pcoupl == 'xyz':
        print('Anisotropic barostat')
        system.addForce(
            mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, True, True, True, nsteps))
    elif pcoupl == 'xy':
        print('Anisotropic barostat only for X and Y')
        system.addForce(
            mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, True, True, False, nsteps))
    elif pcoupl == 'z':
        print('Anisotropic barostat only for Z')
        system.addForce(
            mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, False, False, True, nsteps))
    else:
        raise Exception('Available pressure coupling types: iso, semi-iso, xyz, xy, z')
