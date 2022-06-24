import pytest
import numpy as np
from mstools.topology.geometry import *


def test_grow_particle():
    xyz1 = np.array([0, 0, 0])
    xyz2 = np.array([0, 1, 0])
    print(grow_particle(xyz1, xyz2, 2, np.pi * 0.1))


def test_cluster():
    elements = list(range(10))
    bonds = [(7, 1), (1, 0), (3, 4), (5, 6), (4, 7)]
    matrix = np.zeros((10, 10))
    for i, j in bonds:
        matrix[i][j] = 1
        matrix[j][i] = 1

    assert find_clusters(elements, lambda x, y: matrix[x][y]) == [[0, 1, 3, 4, 7], [2], [5, 6], [8], [9]]
    assert find_clusters_consecutive(elements, lambda x, y: matrix[x][y]) == [[0, 1, 2, 3, 4, 5, 6, 7], [8], [9]]
