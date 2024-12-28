import numpy as np
import openmm
from logging import Logger
from openmm import openmm as mm, app
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
    Run energy minimization on a Simulation until the force on every atom is smaller than tolerance.

    Note that the tolerance here differs from that in OpenMM LocalEnergyMinimizer.
    In OpenMM, the tolerance is set for the root-mean-square of N*3 force components in the whole system.

    For system with constraints, the max force may never reach the tolerance.
    This is because OpenMM uses harmonic restraints instead of rigorous constrints during energy minimization.

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
        if fmax < tolerance or i1_iter >= 5:
            break
        tol_iter /= 4

    state = sim.context.getState(getPositions=True, getEnergy=True)
    if logger:
        logger.info(f'Minimized energy = {state.getPotentialEnergy()}')

    if gro_out:
        GroFile.writeFile(sim.topology, state.getPositions(), state.getPeriodicBoxVectors(), gro_out)


def apply_mc_barostat(system, type, P, T, nstep=100, rigid_molecule=True, logger=None):
    '''
    Add a MonteCarlo barostat to the system and return the added barostat

    Notes
    -----
    It requires the OpenMM with barostat patched https://github.com/openmm/openmm/pull/3797

    Parameters
    ----------
    system : mm.System
    type : str
        The type of barostat.
        Can be one of the following: iso, aniso, xyz-coupled, xy-coupled-z, xyz, xy-coupled, xy, z.
        iso is alias of xyz-coupled; aniso is alias of xyz
    P : float
        Pressure in bar
    T : float
        Temperature in Kelvin
    nstep : int
    rigid_molecule : bool
    logger : Logger

    Returns
    -------
    barostat : mm.Force
        The added barostat force
    '''
    if type in ('xyz-coupled', 'iso'):
        msg = 'Isotropic barostat'
        force = mm.MonteCarloBarostat(P * bar, T * kelvin, nstep, rigid_molecule)
    elif type == 'xy-coupled-z':
        msg = 'Anisotropic barostat with XY coupled'
        force = mm.MonteCarloMembraneBarostat(P * bar, 0 * bar * nm, T * kelvin,
                                              mm.MonteCarloMembraneBarostat.XYIsotropic,
                                              mm.MonteCarloMembraneBarostat.ZFree,
                                              nstep, rigid_molecule)
    elif type in ('xyz', 'aniso'):
        msg = 'Anisotropic barostat'
        force = mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, True, True, True, nstep, rigid_molecule)
    elif type == 'xy-coupled':
        msg = 'Anisotropic barostat with XY coupled and Z fixed'
        force = mm.MonteCarloMembraneBarostat(P * bar, 0 * bar * nm, T * kelvin,
                                              mm.MonteCarloMembraneBarostat.XYIsotropic,
                                              mm.MonteCarloMembraneBarostat.ZFixed,
                                              nstep, rigid_molecule)
    elif type == 'xy':
        msg = 'Anisotropic barostat only for X and Y'
        force = mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, True, True, False, nstep, rigid_molecule)
    elif type == 'z':
        msg = 'Anisotropic barostat only for Z'
        force = mm.MonteCarloAnisotropicBarostat([P * bar] * 3, T * kelvin, False, False, True, nstep, rigid_molecule)
    else:
        raise Exception(
            'Available pressure coupling types: iso, aniso, xyz-coupled, xy-coupled-z, xyz, xy-coupled, xy, z')

    scaled_unit = 'molecules' if rigid_molecule else 'constrained groups'
    msg += f' with {scaled_unit} scaled'

    if logger:
        logger.info(msg)

    system.addForce(force)

    return force


def energy_decomposition(sim: app.Simulation, logger=None):
    '''
    Get the energy contribution of different force groups

    Parameters
    ----------
    sim : app.Simulation
    logger : Logger

    Returns
    -------
    energies : Dict[str, float]
        The name and energy contribution of each force group
    '''
    groups = {}  # {id_group: {force1_name, force2_name, ...}, ...}
    for force in sim.system.getForces():
        i = force.getForceGroup()
        if i not in groups:
            groups[i] = set()
        groups[i].add(force.getName())

    energies = {}
    for i in sorted(groups.keys()):
        energy = sim.context.getState(getEnergy=True, groups={i}).getPotentialEnergy()._value
        if energy != 0 or i != 0:
            name = f'{i}_{"+".join(groups[i])}'
            energies[name] = energy

    if logger:
        for k, v in energies.items():
            logger.info(f'E_{k} = {v}')

    return energies
