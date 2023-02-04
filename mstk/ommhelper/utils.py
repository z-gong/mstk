import openmm
import openmm.openmm as mm
import numpy as np
from openmm import app
from .grofile import GroFile
from .unit import kelvin, bar, nm, ps


class CONST:
    PI = 3.1415926535
    EPS0 = 8.8541878128E-12  # F/m
    ONE_4PI_EPS0 = 138.935455  # qq/r -> kJ/mol


def print_omm_info(logger=None):
    messages = ['OpenMM verison: ' + openmm.version.version,
                'OpenMM library: ' + openmm.version.openmm_library_path,
                'OpenMM platforms: ' + ', '.join([mm.Platform.getPlatform(i).getName()
                                                  for i in range(mm.Platform.getNumPlatforms())]),
                ]
    for msg in messages:
        if logger:
            logger.info(msg)
        else:
            print(msg)


def minimize(sim, tolerance, gro_out=None, logger=None):
    '''
    Run energy minimization on a Simulation until the force on each atom is smaller than tolerance.

    Note that the tolerance here differs from that in OpenMM LocalEnergyMinimizer.
    In OpenMM, the tolerance is set for the root mean square of N*3 force components in the whole system.

    Parameters
    ----------
    sim : app.Simulation
    tolerance : float
    gro_out : str
    logger : Logger
    '''
    if logger:
        logger.info(f'Initial energy = {sim.context.getState(getEnergy=True).getPotentialEnergy()}')
    i1_iter = 0
    tol_iter = tolerance
    while True:
        i1_iter += 1
        sim.minimizeEnergy(tolerance=tol_iter)
        forces = sim.context.getState(getForces=True).getForces(asNumpy=True)._value
        fsq = np.sum(forces * forces, axis=1)
        imax = np.argmax(fsq)
        fmax = fsq[imax] ** 0.5
        if logger:
            logger.info(f'Iter {i1_iter} Max force F_{imax} = {fmax} kJ/mol/nm')
        if fmax < tolerance or i1_iter >= 10:
            break
        tol_iter /= 2

    state = sim.context.getState(getPositions=True, getEnergy=True)
    if logger:
        logger.info(f'Minimized energy = {state.getPotentialEnergy()}')

    if gro_out:
        GroFile.writeFile(sim.topology, state.getPositions(), state.getPeriodicBoxVectors(), gro_out)


def apply_mc_barostat(system, pcoupl, P, T, nsteps=100, logger=None):
    '''
    Add a MonteCarlo barostat to the system and return the index of the added barostat

    Returns
    -------
    idx : index of the added barostat force
    '''
    if pcoupl == 'iso':
        msg = 'Isotropic barostat'
        force = mm.MonteCarloBarostat(P * bar, T * kelvin, nsteps)
    elif pcoupl == 'semi-iso':
        msg = 'Anisotropic barostat with coupled XY'
        force = mm.MonteCarloMembraneBarostat(P * bar, 0 * bar * nm, T * kelvin,
                                              mm.MonteCarloMembraneBarostat.XYIsotropic,
                                              mm.MonteCarloMembraneBarostat.ZFree, nsteps)
    elif pcoupl == 'xyz':
        msg = 'Anisotropic barostat'
        force = mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, True, True, True, nsteps)
    elif pcoupl == 'xy':
        msg = 'Anisotropic barostat only for X and Y'
        force = mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, True, True, False, nsteps)
    elif pcoupl == 'z':
        msg = 'Anisotropic barostat only for Z'
        force = mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, False, False, True, nsteps)
    else:
        raise Exception('Available pressure coupling types: iso, semi-iso, xyz, xy, z')

    if logger:
        logger.info(msg)

    return system.addForce(force)


def energy_decomposition(sim: app.Simulation, logger=None):
    groups = {}
    for force in sim.system.getForces():
        i = force.getForceGroup()
        if i not in groups:
            groups[i] = set()
        groups[i].add(force.getName())

    energies = {}
    for i in sorted(groups.keys()):
        energy = sim.context.getState(getEnergy=True, groups={i}).getPotentialEnergy()._value
        if energy != 0 or i != 0:
            name = '+'.join(groups[i])
            if name in energies:
                name += f'_{i}'
            energies[name] = energy

    if logger:
        for k, v in energies.items():
            logger.info(f'E_{k} = {v}')

    return energies
