import math
import numpy as np


class UnitCell:
    '''
    UnitCell represents the periodic boundary condition of simulation system.

    Both rectangular and triclinic unit cells are supported by this class.
    However, almost all other modules in `mstk` only support rectangular cell.

    Parameters
    ----------
    box : array_like, optional
        box can be array_like of shape (3,), (2,3) or (3,3).
        If its shape is (3,), then it will be treated as the lengths of a rectangular box.
        If its shape is (2,3), then it will be treated as the lengths and angles of a triclinic box.
        The angles are in unit of radian.
        If its shape is (3,3), then it will be treated as the box vectors.
        If not provided, the cell will be initialized with all elements of box vectors equal to zero.
    '''

    def __init__(self, box=None):
        self._vectors = None
        self._lengths = None
        self._angles = None
        self.set_box(box)

    def set_box(self, box):
        '''
        Set the unit cell with lengths or (lengths and angles) or full vectors.

        Parameters
        ----------
        box : array_like or None
            See the argument `box` of constructor for detailed explanations.
        '''
        if box is None:
            if self._vectors is None:
                self._vectors = np.zeros((3, 3), dtype=float)
            else:
                self._vectors.fill(0.0)
            if self._lengths is None:
                self._lengths = np.zeros(3, dtype=float)
            else:
                self._lengths.fill(0.0)
            if self._angles is None:
                self._angles = np.empty(3, dtype=float)
            self._angles.fill(math.pi / 2)
            return

        if not isinstance(box, (tuple, list, np.ndarray)):
            raise ValueError('Invalid argument')
        array = np.array(box, dtype=float)
        if array.shape == (3,):
            self._vectors = np.array([[array[0], 0, 0],
                                      [0, array[1], 0],
                                      [0, 0, array[2]]], dtype=float)
            self._lengths = array
            self._angles = np.array([math.pi / 2] * 3, dtype=float)
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
        '''
        The box vectors of the unit cell.

        The box vectors is a 3*3 matrix [[Vxx, Vxy, Vxz], [Vyx, Vyy, Vyz], [Vzx, Vzy, Vzz]].
        Several criterion will be satisfied:
        Vxy = Vxz = Vyz = 0; Vxx > 2*Vyx; Vxx > 2*Vzx; Vyy > 2*Vzy.

        Returns
        -------
        vectors : array_like
            The shape of box vectors is (3,3)
        '''
        if self._vectors is None:
            self._vectors = self.calc_vectors_from_lengths_angles(self._lengths, self._angles)
        return self._vectors

    def get_size(self):
        '''
        The diagonal values of the box vectors.

        Returns
        -------
        size : array_like
            The shape of the size vector is (3,)

        '''
        return np.array([self.vectors[0][0], self.vectors[1][1], self.vectors[2][2]])

    @property
    def volume(self):
        '''
        The volume of the unit cell.

        Returns
        -------
        volume : float
        '''
        return self.vectors[0][0] * self.vectors[1][1] * self.vectors[2][2]

    @property
    def lengths(self):
        '''
        The lengths of three edges of the unit cell.

        Returns
        -------
        lengths : array_like
            The shape of the length vector is (3,)

        '''
        if self._lengths is None:
            self._lengths, self._angles = self.calc_lengths_angles_from_vectors(self.vectors)
        return self._lengths

    @property
    def angles(self):
        '''
        The angles formed between the three edges of the unit cell, in unit of radian.

        Returns
        -------
        angles: array_like
            The shape of the angle vector is (3,)
        '''

        if self._angles is None:
            self._lengths, self._angles = self.calc_lengths_angles_from_vectors(self.vectors)
        return self._angles

    @property
    def is_rectangular(self):
        '''
        Whether or not this unit cell is rectangular.

        Returns
        -------
        is : bool
        '''
        return all(np.abs(self.angles - math.pi / 2) < 1e-4)

    @staticmethod
    def _reduce_box_vectors(vectors):
        '''
        Transform the orientations so that a0 >= 2*b0, a0 >= 2*c0, b1 >= 2*c1
        '''
        a, b, c = vectors
        if a[0] >= 2 * b[0] and a[0] >= 2 * c[0] and b[1] >= 2 * c[1]:
            return
        c -= b * round(c[1] / b[1])
        c -= a * round(c[0] / a[0])
        b -= a * round(b[0] / a[0])

    @staticmethod
    def calc_vectors_from_lengths_angles(lengths, angles):
        '''
        Convert box lengths and angles to periodic box vectors.

        Parameters
        ----------
        lengths : array_like
            Length vector in shape of (3,)
        angles : array_like
            Angle vector in shape of (3,) in unit o radian

        Returns
        -------
        vectors : array_like
            Box vectors in shape of (3,3)
        '''
        la, lb, lc = lengths
        alpha, beta, gamma = angles

        if all(np.abs(angles - math.pi) / 2 < 1e-4):
            return np.array([[la, 0, 0],
                             [0, lb, 0],
                             [0, 0, lc]], dtype=float)

        a = np.array([la, 0, 0], dtype=float)
        b = np.array([lb * math.cos(gamma), lb * math.sin(gamma), 0], dtype=float)
        cx = lc * math.cos(beta)
        cy = lc * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
        cz = math.sqrt(lc * lc - cx * cx - cy * cy)
        c = np.array([cx, cy, cz], dtype=float)

        # if any element is very close to 0, set it to exactly 0
        for i in range(3):
            if abs(a[i]) < 1e-4:
                a[i] = 0.0
            if abs(b[i]) < 1e-4:
                b[i] = 0.0
            if abs(c[i]) < 1e-4:
                c[i] = 0.0

        # make sure ax > 2*bx, ax > 2*cz, by > 2*bz
        c = c - b * round(c[1] / b[1])
        c = c - a * round(c[0] / a[0])
        b = b - a * round(b[0] / a[0])

        return np.array([a, b, c])

    @staticmethod
    def calc_lengths_angles_from_vectors(vectors):
        '''
        Convert periodic box vectors to lengths and angles.

        Parameters
        ----------
        vectors : array_like
            Box vectors in shape of (3,3)

        Returns
        -------
        lengths : array_like
            Length vector in shape of (3,)
        angles : array_like
            Angle vector in shape of (3,) in unit of radian
        '''
        a, b, c = vectors
        la = np.sqrt(np.dot(a, a))
        lb = np.sqrt(np.dot(b, b))
        lc = np.sqrt(np.dot(c, c))

        if b[1] == 0 and c[0] == 0 and c[1] == 0:
            alpha = beta = gamma = math.pi / 2
        else:
            alpha = math.acos(np.clip(np.dot(b, c) / (lb * lc), -1, 1))
            beta = math.acos(np.clip(np.dot(c, a) / (lc * la), -1, 1))
            gamma = math.acos(np.clip(np.dot(a, b) / (la * lb), -1, 1))

            # if any angle is very close to PI/2, set it to exactly PI/2
            if abs(alpha - math.pi / 2) < 1e-4:
                alpha = math.pi / 2
            if abs(beta - math.pi / 2) < 1e-4:
                beta = math.pi / 2
            if abs(gamma - math.pi / 2) < 1e-4:
                gamma = math.pi / 2

        return np.array([la, lb, lc], dtype=float), np.array([alpha, beta, gamma], dtype=float)
