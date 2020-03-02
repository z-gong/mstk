import simtk.openmm as mm
from simtk import unit
from simtk.unit import elementary_charge as e, nanometer as nm, kilojoule_per_mole as kJ_mol

class CONST:
    PI = 3.1415926535
    EPS0 = 8.8541878128E-12 * unit.farad / unit.meter # vacuum dielectric constant

def slab_correction(system: mm.System):
    '''
    This applies Yeh's long range coulomb correction for slab geometry in z direction
    to eliminate the undesired interactions between periodic slabs
    It's useful for 2-D systems simulated under 3-D periodic condition
    Note that a vacuum space two times larger than slab thickness is required
    for this correction to work correctly
    The box size should not change during the simulation
    '''
    muz = mm.CustomExternalForce('q*z')
    muz.addPerParticleParameter('q')
    nbforce = [f for f in system.getForces() if f.__class__ == mm.NonbondedForce][0]
    qsum = 0
    for i in range(nbforce.getNumParticles()):
        q = nbforce.getParticleParameters(i)[0].value_in_unit(e)
        muz.addParticle(i, [q])
        qsum += q
    if abs(qsum) > 1E-4:
        raise Exception('Slab correction is not valid for charged system')

    box = system.getDefaultPeriodicBoxVectors()
    vol = (box[0][0] * box[1][1] * box[2][2]).value_in_unit(nm ** 3)
    # convert from e^2/nm to kJ/mol  # 138.93545915168772
    _conv = (1 / (4 * CONST.PI * CONST.EPS0) * e**2 / nm / unit.item).value_in_unit(unit.kilojoule_per_mole)
    prefactor = 2 * CONST.PI / vol * _conv
    cvforce = mm.CustomCVForce(f'{prefactor}*muz*muz')
    cvforce.addCollectiveVariable('muz', muz)
    system.addForce(cvforce)

    return cvforce


def spring_self(system: mm.System, positions: [mm.Vec3], particles: [int], kx, ky, kz):
    '''
    Restrain the particles at its original positions
    Note that the original positions will NOT change with box size if there is barostat
    '''
    if system.getNumParticles() != len(positions):
        raise Exception('Length of positions does not equal to number of particles in system')

    if unit.is_quantity(kx):
        kx = kx.value_in_unit(kJ_mol / nm**2)
    if unit.is_quantity(ky):
        ky = ky.value_in_unit(kJ_mol / nm**2)
    if unit.is_quantity(kz):
        kz = kz.value_in_unit(kJ_mol / nm**2)

    force = mm.CustomExternalForce('kx*periodicdistance(x,0,0,x0,0,0)^2+'
                                   'ky*periodicdistance(0,y,0,0,y0,0)^2+'
                                   'kz*periodicdistance(0,0,z,0,0,z0)^2')
    force.addGlobalParameter('kx', kx) # kJ/mol.nm^2
    force.addGlobalParameter('ky', ky)
    force.addGlobalParameter('kz', kz)
    force.addPerParticleParameter('x0')
    force.addPerParticleParameter('y0')
    force.addPerParticleParameter('z0')
    for i in particles:
        force.addParticle(i, list(positions[i]))
    system.addForce(force)

    return force


def wall_power(system: mm.System, particles: [int], direction: str, bound: [float], k, cutoff, power=2):
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
    force = mm.CustomExternalForce(f'k*step({_min_0}-{direction})*rmin^{power}+'
                                   f'k*step({direction}-{_max_0})*rmax^{power};'
                                   f'rmin=({_min_0}-{direction})/{cutoff};'
                                   f'rmax=({direction}-{_max_0})/{cutoff}')
    force.addGlobalParameter('k', k) # kJ/mol
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
    force = mm.CustomExternalForce(f'4*eps*step({_min_0}-{direction})*(rmin^12-rmin^6+0.25)+'
                                   f'4*eps*step({direction}-{_max_0})*(rmax^12-rmax^6+0.25);'
                                   f'rmin={sigma}/({direction}-{_min});'
                                   f'rmax={sigma}/({_max}-{direction})')
    force.addGlobalParameter('eps', epsilon) # kJ/mol
    for i in particles:
        force.addParticle(i, [])
    system.addForce(force)

    return force

def electric_field(system: mm.System, particles: [int], efx, efy, efz):
    '''
    Apply external electric field to particles
    The unit of electric field strength is V/nm
    '''
    if unit.is_quantity(efx):
        efx = efx.value_in_unit(unit.volt / unit.nanometer)
    if unit.is_quantity(efy):
        efy = efy.value_in_unit(unit.volt / unit.nanometer)
    if unit.is_quantity(efz):
        efz = efz.value_in_unit(unit.volt / unit.nanometer)

    # convert from eV/nm to kJ/(mol*nm)  # 96.4853400990037
    _conv = (1 * e * unit.volt / unit.item).value_in_unit(kJ_mol)
    force = mm.CustomExternalForce(f'{_conv}*(efx*q*x+efy*q*y+efz*q*z)')
    force.addGlobalParameter('efx', efx)
    force.addGlobalParameter('efy', efy)
    force.addGlobalParameter('efz', efz)
    force.addPerParticleParameter('q')
    nbforce = [f for f in system.getForces() if f.__class__ == mm.NonbondedForce][0]
    for i in particles:
        q = nbforce.getParticleParameters(i)[0].value_in_unit(e)
        force.addParticle(i, [q])
    system.addForce(force)

    return force
