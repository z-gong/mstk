import math
import numpy as np


class UnitCell():
    def __init__(self, arg):
        self.init_vectors(arg)

    def init_vectors(self, arg):
        if not isinstance(arg, (tuple, list, np.ndarray)) or len(arg) != 3:
            raise ValueError('Invalid argument')
        array = np.array(arg, dtype=np.float32)
        if array.shape == (3, 3):
            self._vectors = array
            self._lengths = None  # calculate lengths and angles only when needed
            self._angles = None  # calculate lengths and angles only when needed
        elif array.shape == (3,):
            self._vectors = np.array([[array[0], 0, 0],
                                      [0, array[1], 0],
                                      [0, 0, array[2]]], dtype=np.float32)
            self._lengths = array
            self._angles = np.array([90, 90, 90], dtype=np.float32)
        else:
            raise ValueError('Invalid argument')

    @property
    def vectors(self):
        return self._vectors

    @vectors.setter
    def vectors(self, value):
        self.init_vectors(value)

    @property
    def box(self):
        return np.array([self.vectors[0][0], self.vectors[1][1], self.vectors[2][2]])

    @property
    def volume(self):
        return self.vectors[0][0] * self.vectors[1][1] * self.vectors[2][2]

    @property
    def lengths(self):
        if self._lengths is None:
            self._lengths, self._angles = self.calc_lengths_angles_from_vectors(self.vectors)
        return self._lengths

    @property
    def angles(self):
        if self._angles is None:
            self._lengths, self._angles = self.calc_lengths_angles_from_vectors(self.vectors)
        return self._angles

    @property
    def is_rectangular(self):
        return all(self.angles == 90)

    @staticmethod
    def calc_vectors_from_lengths_angles(lengths, angles):
        '''
        Convert lengths and angles to periodic box vectors
        Angles should be in degrees
        '''
        la, lb, lc = lengths
        alpha, beta, gamma = angles
        alpha *= math.pi / 180
        beta *= math.pi / 180
        gamma *= math.pi / 180

        a = np.array([la, 0, 0], dtype=np.float32)
        b = np.array([lb * math.cos(gamma), lb * math.sin(gamma), 0], dtype=np.float32)
        cx = lc * math.cos(beta)
        cy = lc * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
        cz = math.sqrt(lc * lc - cx * cx - cy * cy)
        c = np.array([cx, cy, cz], dtype=np.float32)

        # if any element is very close to 0, set it to exactly 0
        for i in range(3):
            if abs(a[i]) < 1e-6:
                a[i] = 0.0
            if abs(b[i]) < 1e-6:
                b[i] = 0.0
            if abs(c[i]) < 1e-6:
                c[i] = 0.0

        # make sure ax > 2*bx, ax > 2*cz, by > 2*bz
        c = c - b * round(c[1] / b[1])
        c = c - a * round(c[0] / a[0])
        b = b - a * round(b[0] / a[0])

        return np.array([a, b, c])

    @staticmethod
    def calc_lengths_angles_from_vectors(vectors):
        '''
        Convert periodic box vectors to lengths and angles
        Angles are returned in degrees
        '''
        a, b, c = vectors[0], vectors[1], vectors[2]
        la = np.sqrt(np.dot(a, a))
        lb = np.sqrt(np.dot(b, b))
        lc = np.sqrt(np.dot(c, c))
        alpha = math.acos(np.clip(np.dot(b, c) / (lb * lc), -1, 1)) * 180 / math.pi
        beta = math.acos(np.clip(np.dot(c, a) / (lc * la), -1, 1)) * 180 / math.pi
        gamma = math.acos(np.clip(np.dot(a, b) / (la * lb), -1, 1)) * 180 / math.pi

        # if any angle is very close to 90 degree, set it to exactly 90 degree
        if abs(alpha - 90) < 1e-6:
            alpha = 90
        if abs(beta - 90) < 1e-6:
            beta = 90
        if abs(gamma - 90) < 1e-6:
            gamma = 90

        return np.array([la, lb, lc], dtype=np.float32), \
               np.array([alpha, beta, gamma], dtype=np.float32)
