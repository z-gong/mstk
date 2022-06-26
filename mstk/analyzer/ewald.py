import math
import numpy as np
from mstk.chem.constant import *


class EwaldSum():
    def __init__(self, charges, positions, box, cutoff=1.2, tolerance=5E-4):
        self.charges = np.array(charges)
        self.positions = np.array(positions)
        self.n_atom = len(self.charges)
        self.box = np.array(box)
        self.volume = box[0] * box[1] * box[2]
        self.reciprocal_box = 2 * PI / self.box
        self.cutoff = cutoff
        self.tolerance = tolerance

        self.alpha, self.kmax_x, self.kmax_y, self.kmax_z = self.calc_ewald_parameters()
        self.eikr = self.calc_eikr()  # exp(ikr)
        self.eikr_conj = np.conj(self.eikr)

    def calc_ewald_parameters(self):
        def get_kmax(width, alpha, tolerance, guess):
            def get_ewald_error(kmax):
                temp = kmax * PI / (width * alpha)
                return 0.05 * kmax * math.sqrt(width * alpha) * math.exp(-temp * temp)

            kmax = guess
            error = get_ewald_error(kmax)
            if error < tolerance:
                while error < tolerance and kmax > 0:
                    kmax -= 1
                    error = get_ewald_error(kmax)
                return kmax + 1
            while error > tolerance:
                kmax += 1
                error = get_ewald_error(kmax)
            return kmax

        lx, ly, lz = self.box
        tol = self.tolerance
        alpha = math.sqrt(-math.log(2.0 * tol)) / self.cutoff
        kmaxx = get_kmax(lx, alpha, tol, 10)
        kmaxy = get_kmax(ly, alpha, tol, 10)
        kmaxz = get_kmax(lz, alpha, tol, 10)
        if kmaxx % 2 == 0: kmaxx += 1
        if kmaxy % 2 == 0: kmaxy += 1
        if kmaxz % 2 == 0: kmaxz += 1

        return alpha, kmaxx, kmaxy, kmaxz

    def calc_eikr(self):
        kmax = max(self.kmax_x, self.kmax_y, self.kmax_z)
        array = np.zeros((kmax, self.n_atom, 3), dtype=complex)
        for n in range(self.n_atom):
            for d in range(3):
                array[0][n][d] = complex(1, 0)
                array[1][n][d] = complex(math.cos(self.positions[n][d] * self.reciprocal_box[d]),
                                         math.sin(self.positions[n][d] * self.reciprocal_box[d]))
                for k in range(2, kmax):
                    array[k][n][d] = array[k - 1][n][d] * array[1][n][d]
        return array

    def distance(self, pos1, pos2, cutoff=None):
        delta = pos2 - pos1
        if cutoff is None:
            return delta, np.sqrt(delta.dot(delta))

        delta = delta - np.floor(delta / self.box + 0.5) * self.box
        if np.any(abs(delta) > cutoff):
            return delta, 0
        r = np.sqrt(delta.dot(delta))
        r = 0 if r > cutoff else r
        return delta, r

    def calc_energy_forces(self, ewald=True):
        if ewald:
            e_short, f_short = self.calc_ewald_short()
            e_long, f_long = self.calc_ewald_long()
            e_self = self.calc_ewald_self()
            return e_short + e_long + e_self, f_short + f_long

        energy = 0
        forces = np.zeros((self.n_atom, 3))
        for i in range(self.n_atom):
            for j in range(i + 1, self.n_atom):
                delta, r = self.distance(self.positions[i], self.positions[j])
                e = ONE_4PI_EPS0 * self.charges[i] * self.charges[j] / r
                energy += e
                for d in range(3):
                    dEdR = e / r / r * delta[d]
                    forces[i][d] -= dEdR
                    forces[j][d] += dEdR
        return energy, forces

    def calc_ewald_short(self):
        energy = 0
        forces = np.zeros((self.n_atom, 3))
        for i in range(self.n_atom):
            for j in range(i + 1, self.n_atom):
                delta, r = self.distance(self.positions[i], self.positions[j], self.cutoff)
                if r > 0:
                    e = ONE_4PI_EPS0 * self.charges[i] * self.charges[j] / r
                    alpha_r = self.alpha * r
                    energy += e * math.erfc(alpha_r)
                    for d in range(3):
                        dEdR = e / r / r * delta[d] * math.erfc(alpha_r) \
                               + e * (-2 / PI_SQRT * math.exp(-alpha_r * alpha_r)) \
                               * self.alpha * (-delta[d] / r)
                        forces[i][d] -= dEdR
                        forces[j][d] += dEdR
        return energy, forces

    def calc_ewald_long(self):
        factor_ewald = -1 / (4 * self.alpha ** 2)
        factor_energy = 1 / (2 * self.volume * VACUUM_PERMITTIVITY) \
                        * ELEMENTARY_CHARGE ** 2 / NANO / 1000 * AVOGADRO
        energy = 0
        forces = np.zeros((self.n_atom, 3))
        qxyz = np.zeros(self.n_atom, dtype=complex)
        # for ix in range(-self.kmax_x + 1, self.kmax_x):
        #     for iy in range(-self.kmax_y + 1, self.kmax_y):
        #         for iz in range(-self.kmax_z + 1, self.kmax_z):
        iy_low = 0
        iz_low = 1
        for ix in range(0, self.kmax_x):
            kx = ix * self.reciprocal_box[0]
            for iy in range(iy_low, self.kmax_y):
                ky = iy * self.reciprocal_box[1]
                for iz in range(iz_low, self.kmax_z):
                    kz = iz * self.reciprocal_box[2]
                    S_k = complex(0, 0)
                    for j in range(self.n_atom):
                        eikr_x = self.eikr[ix][j][0] if ix >= 0 else self.eikr_conj[-ix][j][0]
                        eikr_y = self.eikr[iy][j][1] if iy >= 0 else self.eikr_conj[-iy][j][1]
                        eikr_z = self.eikr[iz][j][2] if iz >= 0 else self.eikr_conj[-iz][j][2]
                        qxyz[j] = self.charges[j] * eikr_x * eikr_y * eikr_z
                        S_k += qxyz[j]
                    k2 = kx * kx + ky * ky + kz * kz
                    ak = math.exp(k2 * factor_ewald) / k2
                    energy += ak * (S_k.real ** 2 + S_k.imag ** 2)
                    for j in range(self.n_atom):
                        # force = -1/(2V*eps0)*ak*d(S(k)S(-k))/dRj
                        # f*k=dS(k)/dRj*S(-k)+S(K)*dS(-k)/dRj
                        f = 2 * (qxyz[j].real * S_k.imag - qxyz[j].imag * S_k.real)
                        forces[j][0] -= ak * f * kx
                        forces[j][1] -= ak * f * ky
                        forces[j][2] -= ak * f * kz
                iz_low = 1 - self.kmax_z
            iy_low = 1 - self.kmax_y
        return energy * factor_energy * 2, forces * factor_energy * 2

    def calc_ewald_self(self):
        return - ONE_4PI_EPS0 * (self.charges * self.charges).sum() * self.alpha / PI_SQRT

    @staticmethod
    def create_test(n_atom=100, box=None, seed=0, **kwargs):
        if box is None:
            box = [3.0, 4.0, 5.0]
        np.random.seed(seed)
        charges = np.random.random(n_atom) - 0.5
        charges[-1] = 0 - charges[:-1].sum()
        positions = np.random.random((n_atom, 3)) * 3
        ewald = EwaldSum(charges, positions, box, **kwargs)
        return ewald

    def calc_energy_forces_with_omm(self, ewald=True, short=True, long=True):
        from openmm import openmm as mm, unit
        system = mm.System()
        system.setDefaultPeriodicBoxVectors(*np.diag(self.box))
        nbforce = mm.NonbondedForce()
        system.addForce(nbforce)

        if ewald:
            nbforce.setNonbondedMethod(mm.NonbondedForce.Ewald)
        else:
            nbforce.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        nbforce.setCutoffDistance(self.cutoff)
        nbforce.setForceGroup(1)
        nbforce.setReciprocalSpaceForceGroup(2)

        for i in range(self.n_atom):
            system.addParticle(0)
            nbforce.addParticle(self.charges[i], 1.0, 0)

        integrator = mm.VerletIntegrator(0.001)
        platform = mm.Platform.getPlatformByName('Reference')
        context = mm.Context(system, integrator, platform)
        context.setPositions(self.positions)
        groups = set()
        if short: groups.add(1)
        if long: groups.add(2)
        state = context.getState(getEnergy=True, getForces=True, groups=groups)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        forces = state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer)
        return energy, forces
