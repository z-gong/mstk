import simtk.openmm as mm
from .util import CONST
from .unit import *


def slab_correction(system: mm.System):
    '''
    This applies Yeh's long range coulomb correction for slab geometry in z direction
    to eliminate the undesired interactions between periodic slabs
    It's useful for 2-D systems simulated under 3-D periodic condition
    For this correction to work correctly:
    * A vacuum space two times larger than slab thickness is required
    * All particles should never diffuse across the z boundaries
    * The box size should not change during the simulation
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
        raise Exception('Slab correction is not valid for charged system')

    box = system.getDefaultPeriodicBoxVectors()
    vol = (box[0][0] * box[1][1] * box[2][2]).value_in_unit(nm ** 3)
    # convert from e^2/nm to kJ/mol  # 138.93545915168772
    _eps0 = CONST.EPS0 * unit.farad / unit.meter  # vacuum dielectric constant
    _conv = (1 / (4 * CONST.PI * _eps0) * qe ** 2 / nm / unit.item).value_in_unit(kJ_mol)
    prefactor = 2 * CONST.PI / vol * _conv
    cvforce = mm.CustomCVForce(f'{prefactor}*muz*muz')
    cvforce.addCollectiveVariable('muz', muz)
    system.addForce(cvforce)

    return cvforce


def spring_self(system: mm.System, positions: [mm.Vec3], particles: [int], strength):
    '''
    Restrain the particles at its original positions
    Note that the original positions will NOT change with box size if there is barostat
    '''
    if system.getNumParticles() != len(positions):
        raise Exception('Length of positions does not equal to number of particles in system')

    if unit.is_quantity(strength):
        kx, ky, kz = strength.value_in_unit(kJ_mol / nm ** 2)
    else:
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


def wall_power(system: mm.System, particles: [int], direction: str, bound: [float],
               k, cutoff, power=2):
    '''
    Set a wall for particles so that they cannot cross it
    Note that periodic box condition is not considered,
    so you need to make sure particles will not move to other cells during the simulation

    The energy equal to k when particle is located at the lower or higher bound
    and equal to zero when particle is located between [lower bound + cutoff, higher bound - cutoff]
    '''
    if direction not in ['x', 'y', 'z']:
        raise Exception('direction can only be x, y or z')
    _min, _max = bound
    if unit.is_quantity(_min):
        _min = _min.value_in_unit(nm)
    if unit.is_quantity(_max):
        _max = _max.value_in_unit(nm)
    if unit.is_quantity(k):
        k = k.value_in_unit(kJ_mol)
    if unit.is_quantity(cutoff):
        cutoff = cutoff.value_in_unit(nm)

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


def wall_lj126(system: mm.System, particles: [int], direction: str, bound: [float], epsilon, sigma):
    '''
    Set a wall for particles so that they cannot cross it
    Note that periodic box condition is not considered,
    so you need to make sure particles will not move to other cells during the simulation

    The energy is infinite when particle is located at the lower or higher bound
    and equal to epsilon when particle is located at lower bound + sigma or higher bound - sigma
    and equal to zero when particle is located between [lower bound + sigma * 2^(1/6), higher bound - sigma * 2^(1/6)]
    '''
    if direction not in ['x', 'y', 'z']:
        raise Exception('direction can only be x, y or z')
    _min, _max = bound
    if unit.is_quantity(_min):
        _min = _min.value_in_unit(nm)
    if unit.is_quantity(_max):
        _max = _max.value_in_unit(nm)
    if unit.is_quantity(epsilon):
        epsilon = epsilon.value_in_unit(kJ_mol)
    if unit.is_quantity(sigma):
        sigma = sigma.value_in_unit(nm)

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


def electric_field(system: mm.System, particles: [int], strength):
    '''
    Apply external electric field to particles
    The unit of electric field strength is V/nm
    '''
    if unit.is_quantity(strength):
        efx, efy, efz = strength.value_in_unit(unit.volt / unit.nanometer)
    else:
        efx, efy, efz = strength

    # convert from eV/nm to kJ/(mol*nm)  # 96.4853400990037
    _conv = (1 * qe * unit.volt / unit.item).value_in_unit(kJ_mol)
    force = mm.CustomExternalForce(f'{_conv}*({efx}*q*x+{efy}*q*y+{efz}*q*z)')
    force.addPerParticleParameter('q')
    nbforce = [f for f in system.getForces() if f.__class__ == mm.NonbondedForce][0]
    for i in particles:
        q = nbforce.getParticleParameters(i)[0].value_in_unit(qe)
        force.addParticle(i, [q])
    system.addForce(force)

    return force


def CLPolCoulTT(system: mm.System, donors: [int], b: float = 45.0):
    '''
    Apply Tang-Toennies damping between H-bond hydrogen atoms and Drude dipoles
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


def restrain_particle_number(system: mm.System, particles: [int], direction: str, bound,
                             sigma, target, k):
    '''
    Restrain the number of particles in a region
    The region is defined by direction (x, y or z) and bound (lower and upper)
    Each particle is consider as a Gaussian distribution with standard deviation equal to sigma
    The number of particles is restrained to the target value
    using a harmonic function with force constant k
    '''
    if direction not in ['x', 'y', 'z']:
        raise Exception('direction can only be x, y or z')
    _min, _max = bound
    if unit.is_quantity(_min):
        _min = _min.value_in_unit(nm)
    if unit.is_quantity(_max):
        _max = _max.value_in_unit(nm)
    if unit.is_quantity(sigma):
        sigma = sigma.value_in_unit(nm)
    if unit.is_quantity(k):
        k = k.value_in_unit(kJ_mol)

    if _min is not None:
        str_min = f'erf(({_min}-{direction})/{2 ** 0.5 * sigma})'
    else:
        str_min = '-1'

    if _max is not None:
        str_max = f'erf(({_max}-{direction})/{2 ** 0.5 * sigma})'
    else:
        str_max = '1'

    nforce = mm.CustomExternalForce(f'0.5*({str_max}-{str_min})')
    for i in particles:
        nforce.addParticle(i, [])

    cvforce = mm.CustomCVForce(f'0.5*{k}*(number-{target})^2')
    cvforce.addCollectiveVariable('number', nforce)
    system.addForce(cvforce)

    return cvforce
