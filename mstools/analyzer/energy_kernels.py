import numpy as np
from ..chem import constant


def ew_dot(vec1, vec2):
    '''
    Element-wise dot product
    '''
    return np.sum(vec1 * vec2, axis=1)


class HarmonicBondKernel():
    '''
    E = k (r-r0)^2

    Parameters
    ----------
    positions : list of list of float
    indexes : list of list of int
    parameters : list of list of float
    '''

    def __init__(self, positions, indexes, parameters):
        self.positions = np.array(positions, dtype=np.float64)
        self.indexes = np.array(indexes, dtype=int)
        self.a1 = self.indexes[:, 0]
        self.a2 = self.indexes[:, 1]
        self.parameters = np.array(parameters, dtype=np.float64)
        self.r0 = self.parameters[:, 0]
        self.k = self.parameters[:, 1]

    def evaluate(self):
        delta = self.positions[self.a2] - self.positions[self.a1]
        rsq = ew_dot(delta, delta)
        r = np.sqrt(rsq)
        energy = self.k * (r - self.r0) ** 2

        forces_a1 = (2 * self.k * (r - self.r0) / r)[:, np.newaxis] * delta

        forces = np.zeros(self.positions.shape, dtype=np.float64)
        # TODO vectorize this loop
        for i in range(self.positions.shape[0]):
            forces[i] += forces_a1[self.a1 == i].sum(axis=0)
            forces[i] -= forces_a1[self.a2 == i].sum(axis=0)

        return r, energy, forces


class HarmonicAngleKernel():
    '''
    E = k (theta-theta0)^2

    '''

    def __init__(self, positions, indexes, parameters):
        self.positions = np.array(positions, dtype=np.float64)
        self.indexes = np.array(indexes, dtype=int)
        self.a1 = self.indexes[:, 0]
        self.a2 = self.indexes[:, 1]
        self.a3 = self.indexes[:, 2]
        self.parameters = np.array(parameters, dtype=np.float64)
        self.theta0 = self.parameters[:, 0]
        self.k = self.parameters[:, 1]

    def evaluate(self):
        vec1 = self.positions[self.a1] - self.positions[self.a2]
        vec2 = self.positions[self.a3] - self.positions[self.a2]
        r1 = np.sqrt(ew_dot(vec1, vec1))
        r2 = np.sqrt(ew_dot(vec2, vec2))
        cos = ew_dot(vec1, vec2) / r1 / r2
        np.clip(cos, -1, 1, out=cos)
        theta = np.arccos(cos)
        energy = self.k * (theta - self.theta0) ** 2

        sin = np.sqrt(1 - cos * cos)
        sin[sin < 1E-4] = 1E-4

        factor = 2 * self.k * (theta - self.theta0) / sin
        c11 = -factor * cos / r1 / r1
        c12 = factor / r1 / r2
        c31 = -factor * cos / r2 / r2

        forces_a1 = c11[:, np.newaxis] * vec1 + c12[:, np.newaxis] * vec2
        forces_a3 = c31[:, np.newaxis] * vec2 + c12[:, np.newaxis] * vec1

        forces = np.zeros(self.positions.shape, dtype=np.float64)
        # TODO vectorize this loop
        for i in range(self.positions.shape[0]):
            forces[i] += forces_a1[self.a1 == i].sum(axis=0)
            forces[i] += forces_a3[self.a3 == i].sum(axis=0)
            forces[i] -= forces_a1[self.a2 == i].sum(axis=0) + forces_a3[self.a2 == i].sum(axis=0)

        return theta, energy, forces


class OplsTorsionKernel():
    '''
    E = k1 (1+cos(phi)) + k2 (1-cos(2 phi)) + k3 (1+cos(3 phi)) + k4 (1-cos(4 phi))

    '''

    def __init__(self, positions, indexes, parameters):
        self.positions = np.array(positions, dtype=np.float64)
        self.indexes = np.array(indexes, dtype=int)
        self.a1 = self.indexes[:, 0]
        self.a2 = self.indexes[:, 1]
        self.a3 = self.indexes[:, 2]
        self.a4 = self.indexes[:, 3]
        self.parameters = np.array(parameters, dtype=np.float64)
        self.k1 = self.parameters[:, 0]
        self.k2 = self.parameters[:, 1]
        self.k3 = self.parameters[:, 2]
        self.k4 = self.parameters[:, 3]

    def evaluate(self):
        vec1 = self.positions[self.a1] - self.positions[self.a2]
        vec2 = self.positions[self.a3] - self.positions[self.a2]
        vec3 = self.positions[self.a3] - self.positions[self.a4]
        n1 = np.cross(vec1, vec2)
        n2 = np.cross(vec2, vec3)
        rsq_n1 = ew_dot(n1, n1)
        rsq_n2 = ew_dot(n2, n2)
        cos = ew_dot(n1, n2) / np.sqrt(rsq_n1 * rsq_n2)
        np.clip(cos, -1, 1, out=cos)
        sign = np.ones(vec1.shape[0])
        sign[ew_dot(vec1, n2) < 0] = -1
        phi = sign * np.arccos(cos)

        c0 = self.k1 + 2 * self.k2 + self.k3
        c1 = self.k1 - 3 * self.k3
        c2 = -2 * self.k2 + 8 * self.k4
        c3 = 4 * self.k3
        c4 = -8 * self.k4
        cos2 = cos * cos
        cos3 = cos2 * cos
        cos4 = cos3 * cos
        energy = c0 + c1 * cos + c2 * cos2 + c3 * cos3 + c4 * cos4

        sin = np.sin(phi)
        ### (-dE / d phi)
        factor = (c1 * sin
                  + 2 * c2 * cos * sin
                  + 3 * c3 * cos2 * sin
                  + 4 * c4 * cos3 * sin)

        ### (d phi / d ri) copied from OpenMM
        rsq2 = ew_dot(vec2, vec2)
        r2 = np.sqrt(rsq2)

        factor1 = factor * r2 / rsq_n1
        factor4 = -factor * r2 / rsq_n2

        forces_a1 = factor1[:, np.newaxis] * n1
        forces_a4 = factor4[:, np.newaxis] * n2

        factor2 = ew_dot(vec1, vec2) / rsq2
        factor3 = ew_dot(vec3, vec2) / rsq2

        s = factor2[:, np.newaxis] * forces_a1 - factor3[:, np.newaxis] * forces_a4

        forces_a2 = -forces_a1 + s
        forces_a3 = -forces_a4 - s
        ###

        forces = np.zeros(self.positions.shape, dtype=np.float64)
        # TODO vectorize this loop
        for i in range(self.positions.shape[0]):
            forces[i] += forces_a1[self.a1 == i].sum(axis=0)
            forces[i] += forces_a2[self.a2 == i].sum(axis=0)
            forces[i] += forces_a3[self.a3 == i].sum(axis=0)
            forces[i] += forces_a4[self.a4 == i].sum(axis=0)

        return phi, energy, forces


class HarmonicTorsionKernel():
    '''
    E = k (phi-phi0)^2

    '''

    def __init__(self, positions, indexes, parameters):
        self.positions = np.array(positions, dtype=np.float64)
        self.indexes = np.array(indexes, dtype=int)
        self.a1 = self.indexes[:, 0]
        self.a2 = self.indexes[:, 1]
        self.a3 = self.indexes[:, 2]
        self.a4 = self.indexes[:, 3]
        self.parameters = np.array(parameters, dtype=np.float64)
        self.phi0 = self.parameters[:, 0]
        self.k = self.parameters[:, 1]

    def evaluate(self):
        vec1 = self.positions[self.a1] - self.positions[self.a2]
        vec2 = self.positions[self.a3] - self.positions[self.a2]
        vec3 = self.positions[self.a3] - self.positions[self.a4]
        n1 = np.cross(vec1, vec2)
        n2 = np.cross(vec2, vec3)
        rsq_n1 = ew_dot(n1, n1)
        rsq_n2 = ew_dot(n2, n2)
        cos = ew_dot(n1, n2) / np.sqrt(rsq_n1 * rsq_n2)
        np.clip(cos, -1, 1, out=cos)
        sign = np.ones(vec1.shape[0])
        sign[ew_dot(vec1, n2) < 0] = -1
        phi = sign * np.arccos(cos)

        d_phi = phi - self.phi0
        d_phi_abs = np.abs(d_phi)
        d_phi_2pi = 2 * np.pi - d_phi_abs
        d_phi_min = np.min([d_phi_abs, d_phi_2pi], axis=0)
        energy = self.k * d_phi_min ** 2

        ### (-dE / d phi)
        factor = - 2 * self.k * d_phi_min

        sign.fill(1)
        sign[d_phi_2pi < d_phi_abs] = -1
        factor *= sign

        sign.fill(1)
        sign[d_phi < 0] = -1
        factor *= sign

        ### (d phi / d ri) copied from OpenMM
        rsq2 = ew_dot(vec2, vec2)
        r2 = np.sqrt(rsq2)

        factor1 = factor * r2 / rsq_n1
        factor4 = -factor * r2 / rsq_n2

        forces_a1 = factor1[:, np.newaxis] * n1
        forces_a4 = factor4[:, np.newaxis] * n2

        factor2 = ew_dot(vec1, vec2) / rsq2
        factor3 = ew_dot(vec3, vec2) / rsq2

        s = factor2[:, np.newaxis] * forces_a1 - factor3[:, np.newaxis] * forces_a4

        forces_a2 = -forces_a1 + s
        forces_a3 = -forces_a4 - s
        ###

        forces = np.zeros(self.positions.shape, dtype=np.float64)
        # TODO vectorize this loop
        for i in range(self.positions.shape[0]):
            forces[i] += forces_a1[self.a1 == i].sum(axis=0)
            forces[i] += forces_a2[self.a2 == i].sum(axis=0)
            forces[i] += forces_a3[self.a3 == i].sum(axis=0)
            forces[i] += forces_a4[self.a4 == i].sum(axis=0)

        return phi, energy, forces


class ConstrainedTorsionKernel():
    '''
    E = k (1-cos(phi-phi0))

    '''

    def __init__(self, positions, indexes, parameters):
        self.positions = np.array(positions, dtype=np.float64)
        self.indexes = np.array(indexes, dtype=int)
        self.a1 = self.indexes[:, 0]
        self.a2 = self.indexes[:, 1]
        self.a3 = self.indexes[:, 2]
        self.a4 = self.indexes[:, 3]
        self.parameters = np.array(parameters, dtype=np.float64)
        self.phi0 = self.parameters[:, 0]
        self.k = self.parameters[:, 1]

    def evaluate(self):
        vec1 = self.positions[self.a1] - self.positions[self.a2]
        vec2 = self.positions[self.a3] - self.positions[self.a2]
        vec3 = self.positions[self.a3] - self.positions[self.a4]
        n1 = np.cross(vec1, vec2)
        n2 = np.cross(vec2, vec3)
        rsq_n1 = ew_dot(n1, n1)
        rsq_n2 = ew_dot(n2, n2)
        cos = ew_dot(n1, n2) / np.sqrt(rsq_n1 * rsq_n2)
        np.clip(cos, -1, 1, out=cos)
        sign = np.ones(vec1.shape[0])
        sign[ew_dot(vec1, n2) < 0] = -1
        phi = sign * np.arccos(cos)

        d_phi = phi - self.phi0
        energy = self.k * (1 - np.cos(d_phi))

        ### (-dE / d phi)
        factor = -self.k * np.sin(d_phi)

        ### (d phi / d ri) copied from OpenMM
        rsq2 = ew_dot(vec2, vec2)
        r2 = np.sqrt(rsq2)

        factor1 = factor * r2 / rsq_n1
        factor4 = -factor * r2 / rsq_n2

        forces_a1 = factor1[:, np.newaxis] * n1
        forces_a4 = factor4[:, np.newaxis] * n2

        factor2 = ew_dot(vec1, vec2) / rsq2
        factor3 = ew_dot(vec3, vec2) / rsq2

        s = factor2[:, np.newaxis] * forces_a1 - factor3[:, np.newaxis] * forces_a4

        forces_a2 = -forces_a1 + s
        forces_a3 = -forces_a4 - s
        ###

        forces = np.zeros(self.positions.shape, dtype=np.float64)
        # TODO vectorize this loop
        for i in range(self.positions.shape[0]):
            forces[i] += forces_a1[self.a1 == i].sum(axis=0)
            forces[i] += forces_a2[self.a2 == i].sum(axis=0)
            forces[i] += forces_a3[self.a3 == i].sum(axis=0)
            forces[i] += forces_a4[self.a4 == i].sum(axis=0)

        return phi, energy, forces


class NonbondedKernel():
    '''
    E = 4 * eps*((sig/r)^12 - (sig/r)^6) + 138.935455 * qq/r
    '''

    def __init__(self, positions, indexes, parameters):
        self.positions = np.array(positions, dtype=np.float64)
        self.indexes = np.array(indexes, dtype=int)
        self.a1 = self.indexes[:, 0]
        self.a2 = self.indexes[:, 1]
        self.parameters = np.array(parameters, dtype=np.float64)
        self.c12 = 4 * self.parameters[:, 0] * self.parameters[:, 1] ** 12
        self.c6 = 4 * self.parameters[:, 0] * self.parameters[:, 1] ** 6
        self.qq = self.parameters[:, 2]
        self.qqconv = constant.ONE_4PI_EPS0

    def evaluate(self):
        delta = self.positions[self.a2] - self.positions[self.a1]
        rsq = ew_dot(delta, delta)
        r = np.sqrt(rsq)
        r6 = r ** 6
        r12 = r6 * r6
        energy = self.c12 / r12 - self.c6 / r6 + self.qqconv * self.qq / r

        forces_a1 = -(12 * self.c12 / r12 / rsq - 6 * self.c6 / r6 / rsq
                      + self.qqconv * self.qq / rsq / r) \
            [:, np.newaxis] * delta

        forces = np.zeros(self.positions.shape, dtype=np.float64)
        # TODO vectorize this loop
        for i in range(self.positions.shape[0]):
            forces[i] += forces_a1[self.a1 == i].sum(axis=0)
            forces[i] -= forces_a1[self.a2 == i].sum(axis=0)

        return r, energy, forces
