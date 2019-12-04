import simtk.openmm as mm
from simtk.unit import elementary_charge as e, nanometer as nm

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
