import math
import numpy as np


class UnitCell():
    def __init__(self, arg):
        self.set_box(arg)

    def set_box(self, arg):
        if not isinstance(arg, (tuple, list, np.ndarray)):
            raise ValueError('Invalid argument')
        array = np.array(arg, dtype=np.float32)
        if array.shape == (3,):
            self._vectors = np.array([[array[0], 0, 0],
                                      [0, array[1], 0],
                                      [0, 0, array[2]]], dtype=np.float32)
            self._lengths = array
            self._angles = np.array([90, 90, 90], dtype=np.float32)
        elif array.shape == (3, 3):
            if array[0][1] != 0 or array[0][2] != 0 or array[1][2] != 0:
                raise ValueError('Invalid value for box vectors')
            self._reduce_box_vectors(array)
            self._vectors = array
            self._lengths = None  # calculate lengths and angles when required
            self._angles = None  # calculate lengths and angles when required
        elif array.shape == (2, 3):
            self._vectors = None  # calculate vectors when required
            self._lengths = array[0]
            self._angles = array[1]

        else:
            raise ValueError('Invalid argument')

    @property
    def vectors(self):
        if self._vectors is None:
            self._vectors = self.calc_vectors_from_lengths_angles(self._lengths, self._angles)
        return self._vectors

    @property
    def size(self):
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
    def _reduce_box_vectors(vectors):
        """
        Reduces the representation of the PBC
        """
        a, b, c = vectors
        c -= b * round(c[1] / b[1])
        c -= a * round(c[0] / a[0])
        b -= a * round(b[0] / a[0])

    @staticmethod
    def calc_vectors_from_lengths_angles(lengths, angles):
        '''
        Convert lengths and angles to periodic box vectors
        Angles should be in degrees
        '''
        la, lb, lc = lengths
        alpha, beta, gamma = angles

        if alpha == 90 and beta == 90 and gamma == 90:
            return np.array([[la, 0, 0],
                             [0, lb, 0],
                             [0, 0, lc]], dtype=np.float32)

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
        a, b, c = vectors
        la = np.sqrt(np.dot(a, a))
        lb = np.sqrt(np.dot(b, b))
        lc = np.sqrt(np.dot(c, c))

        if b[1] == 0 and c[0] == 0 and c[1] == 0:
            alpha = beta = gamma = 90
        else:
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
