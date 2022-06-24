import openmm.openmm as mm
from .utils import CONST
from .unit import *


def slab_correction(system):
    '''
    Apply Yeh's long range coulomb correction for slab geometry in z direction
    to eliminate the undesired interactions between periodic slabs.

    It's useful for 2-D systems simulated under 3-D periodic condition.
    For this correction to work correctly:

    * A vacuum space two times larger than slab thickness is required.
    * All particles should never diffuse across the z boundaries.
    * The box size should not change during the simulation.

    Parameters
    ----------
    system : mm.System
        The OpenMM system to be simulated

    Returns
    -------
    force : mm.CustomCVForce
    '''
    muz = mm.CustomExternalForce('q*z')
    muz.addPerParticleParameter('q')
    nbforce = [f for f in system.getForces() if f.__class__ == mm.NonbondedForce][0]
    qsum = 0
    for i in range(nbforce.getNumParticles()):
        q = nbforce.getParticleParameters(i)[0].value_in_unit(qe)
        muz.addParticle(i, [q])
        qsum += q
    if abs(qsum) > 1E-4:
        raise Exception('Slab correction is not valid for non-neutral system')

    box = system.getDefaultPeriodicBoxVectors()
    vol = (box[0][0] * box[1][1] * box[2][2]).value_in_unit(nm ** 3)
    # convert from e^2/nm to kJ/mol  # 138.93545915168772
    _eps0 = CONST.EPS0 * farad / meter  # vacuum dielectric constant
    _conv = (1 / (4 * CONST.PI * _eps0) * qe ** 2 / nm / item).value_in_unit(kJ_mol)
    prefactor = 2 * CONST.PI / vol * _conv
    cvforce = mm.CustomCVForce(f'{prefactor}*muz*muz')
    cvforce.addCollectiveVariable('muz', muz)
    system.addForce(cvforce)

    return cvforce


def spring_self(system, positions, particles, strength):
    '''
    Restrain the selected particles at their original positions.

    Note that the original positions will NOT change with box size if there is a barostat.

    Parameters
    ----------
    system : mm.System
        The OpenMM system to be simulated
    positions : array_like
        The positions of all particles in the system in unit of nm
    particles : list of int
        The indexes of particles to be restrained
    strength : list of float
        The strength of harmonic restraint in x, y and z directions.
        Three elements should be provided in unit of kJ/mol/nm^2

    Returns
    -------
    force : mm.CustomExternalForce
    '''
    if system.getNumParticles() != len(positions):
        raise Exception('Length of positions does not equal to number of particles in system')

    kx, ky, kz = strength
    force = mm.CustomExternalForce(f'{kx}*periodicdistance(x,0,0,x0,0,0)^2+'
                                   f'{ky}*periodicdistance(0,y,0,0,y0,0)^2+'
                                   f'{kz}*periodicdistance(0,0,z,0,0,z0)^2')
    force.addPerParticleParameter('x0')
    force.addPerParticleParameter('y0')
    force.addPerParticleParameter('z0')
    for i in particles:
        force.addParticle(i, list(positions[i]))
    system.addForce(force)

    return force


def wall_power(system, particles, direction, bound, k, cutoff, power=2):
    '''
    Add a power wall for selected particles so that they cannot cross it.

    Note that periodic box condition is not considered,
    so you need to make sure particles will not move to other cells during the simulation.

    The energy equal to k when particle is located at the lower or higher bound,
    and equal to zero when particle is located between [lower bound + cutoff, higher bound - cutoff].

    Parameters
    ----------
    system : mm.System
    particles : list of int
    direction : ['x', 'y', 'z']
    bound : list of float
    k : float
    cutoff : float
    power : int, optional

    Returns
    -------
    force : mm.CustomExternalForce

    '''
    if direction not in ['x', 'y', 'z']:
        raise Exception('direction can only be x, y or z')
    _min, _max = bound
    _min_0 = _min + cutoff
    _max_0 = _max - cutoff
    force = mm.CustomExternalForce(f'{k}*step({_min_0}-{direction})*rmin^{power}+'
                                   f'{k}*step({direction}-{_max_0})*rmax^{power};'
                                   f'rmin=({_min_0}-{direction})/{cutoff};'
                                   f'rmax=({direction}-{_max_0})/{cutoff}')
    for i in particles:
        force.addParticle(i, [])
    system.addForce(force)

    return force


def wall_lj126(system, particles, direction, bound, epsilon, sigma):
    '''
    Add a LJ-12-6 wall for selected particles so that they cannot cross it.

    Note that periodic box condition is not considered,
    so you need to make sure particles will not move to other cells during the simulation.

    The energy is infinite when particle is located at the lower or higher bound,
    and equal to epsilon when particle is located at lower bound + sigma or higher bound - sigma,
    and equal to zero when particle is located between [lower bound + sigma * 2^(1/6), higher bound - sigma * 2^(1/6)].

    Parameters
    ----------
    system : mm.System
    particles : list of int
    direction : str
    bound : list of float
    epsilon : float
    sigma : float

    Returns
    -------
    force : mm.CustomExternalForce

    '''
    if direction not in ['x', 'y', 'z']:
        raise Exception('direction can only be x, y or z')
    _min, _max = bound
    _min_0 = _min + sigma * 2 ** (1 / 6)
    _max_0 = _max - sigma * 2 ** (1 / 6)
    force = mm.CustomExternalForce(f'4*{epsilon}*step({_min_0}-{direction})*(rmin^12-rmin^6+0.25)+'
                                   f'4*{epsilon}*step({direction}-{_max_0})*(rmax^12-rmax^6+0.25);'
                                   f'rmin={sigma}/({direction}-{_min});'
                                   f'rmax={sigma}/({_max}-{direction})')
    for i in particles:
        force.addParticle(i, [])
    system.addForce(force)

    return force


def electric_field(system, particles, strength):
    '''
    Apply external electric field to selected particles in a system.

    The unit of electric field strength is V/nm.

    Parameters
    ----------
    system : mm.System
    particles : list of int
    strength : list of float
        Strength of electric field in x, y and z directions, in unit of volt/nm

    Returns
    -------
    force : mm.CustomExternalForce

    '''
    # convert from eV/nm to kJ/(mol*nm)  # 96.4853400990037
    _conv = (1 * qe * volt / item).value_in_unit(kJ_mol)
    efx, efy, efz = strength
    force = mm.CustomExternalForce(f'{_conv}*({efx}*q*x+{efy}*q*y+{efz}*q*z)')
    force.addPerParticleParameter('q')
    nbforce = next(f for f in system.getForces() if type(f) is mm.NonbondedForce)
    for i in particles:
        q = nbforce.getParticleParameters(i)[0].value_in_unit(qe)
        force.addParticle(i, [q])
    system.addForce(force)

    return force


def CLPolCoulTT(system, donors, b=45.0):
    '''
    Apply Tang-Toennies damping for the Coulomb interactions between selected H-bond donors and Drude dipoles.

    Parameters
    ----------
    system : mm.System
    donors : list of int
        Indexes of particles served as H-bond donors
    b : float, optional
        b in unit of /nm

    Returns
    -------
    force : mm.CustomNonbondedForce

    '''
    nbforce: mm.NonbondedForce = next(f for f in system.getForces() if type(f) == mm.NonbondedForce)
    dforce: mm.DrudeForce = next(f for f in system.getForces() if type(f) == mm.DrudeForce)

    drude_pairs = {}  # {parent: drude}
    dipole_set = set()
    for i in range(dforce.getNumParticles()):
        drude, parent, p2, p3, p4, q, alpha, aniso12, aniso34 = dforce.getParticleParameters(i)
        drude_pairs[parent] = drude
        dipole_set.add(parent)
        dipole_set.add(drude)

    ttforce = mm.CustomNonbondedForce('-%.6f*q1*q2/r*beta*gamma;'
                                      'beta=exp(-br);'
                                      'gamma=1+br+br*br/2+br2*br/6+br2*br2/24;'
                                      'br2=br*br;'
                                      'br=%.6f*r' % (CONST.ONE_4PI_EPS0, b))
    ttforce.addPerParticleParameter('q')
    for i in range(system.getNumParticles()):
        if i in drude_pairs:
            q, _, _ = nbforce.getParticleParameters(drude_pairs[i])
            q = -q
        else:
            q, _, _ = nbforce.getParticleParameters(i)
        q = q.value_in_unit(qe)
        ttforce.addParticle([q])
    ttforce.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    ttforce.setCutoffDistance(1.2 * nm)
    ttforce.addInteractionGroup(set(donors), dipole_set)
    ttforce.setForceGroup(9)
    # map all the exclusions from NonbondedForce
    for i in range(nbforce.getNumExceptions()):
        ii, jj, _, _, _ = nbforce.getExceptionParameters(i)
        ttforce.addExclusion(ii, jj)
    system.addForce(ttforce)

    return ttforce


def restrain_particle_number(system, particles, direction, bound,
                             sigma, target, k, weights=None):
    '''
    Restrain the number of selected particles in a region.

    The region is defined by direction (x, y or z) and bound (lower and upper).
    Each particle is consider as a Gaussian distribution with standard deviation equal to sigma.
    The number of particles is restrained to the target value
    using a harmonic function with force constant k.

    Parameters
    ----------
    system : mm.System
    particles : list of int
    direction : ['x', 'y', 'z']
    bound : list of float
    sigma : float
        Variance of the particle Gaussian in unit of nm
    target : float
    k : float
        Strength of the harmonic restraint in unit of kJ/mol
    weights : list of float, optional

    Returns
    -------
    force : mm.CustomCVForce

    '''
    if direction not in ['x', 'y', 'z']:
        raise Exception('direction can only be x, y or z')
    _min, _max = bound
    if weights is None:
        weights = [1.0] * len(particles)
    if len(weights) != len(particles):
        raise Exception('particles and weights should have the same length')

    if _min is not None:
        str_min = f'erf(({_min}-{direction})/{2 ** 0.5 * sigma})'
    else:
        str_min = '-1'

    if _max is not None:
        str_max = f'erf(({_max}-{direction})/{2 ** 0.5 * sigma})'
    else:
        str_max = '1'

    nforce = mm.CustomExternalForce(f'0.5*({str_max}-{str_min})*weight')
    nforce.addPerParticleParameter('weight')
    for i, w in zip(particles, weights):
        nforce.addParticle(i, [w])

    cvforce = mm.CustomCVForce(f'0.5*{k}*(number-{target})^2')
    cvforce.addCollectiveVariable('number', nforce)
    system.addForce(cvforce)

    return cvforce


def point_wall_power(system, particles, reference, k, cutoff, power=2):
    '''
    Add a spherical wall centered at a reference point so that selected particles are repelled from this point.

    The energy equal to k when particle is located at the reference point,
    and equal to zero when particle stays beyond the cutoff.

    Parameters
    ----------
    system : mm.System
    particles : list of int
    reference : list of float
    k : float
    cutoff : float
    power : int, optional

    Returns
    -------
    force : mm.CustomExternalForce
    '''
    x0, y0, z0 = reference
    force = mm.CustomExternalForce(f'step({cutoff}-r)*{k}*(1-r/{cutoff})^{power};'
                                   f'r=periodicdistance(x,y,z,{x0},{y0},{z0});')
    for i in particles:
        force.addParticle(i, [])
    force.setName('PointWall')
    system.addForce(force)

    return force
