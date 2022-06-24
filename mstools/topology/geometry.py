import math
import numpy as np


def transform_coor(unit_vectors, xyz):
    '''
    Transform the coordinates in new system back to its coordinates in old system

    Parameters
    ----------
    unit_vectors : list of array
        unit_vectors of the new system in old system
    xyz : array
        coordinates in new system

    Returns
    -------
    xyz_old : array
        coordinates in old system
    '''
    R = np.array(unit_vectors).T
    b = np.matmul(R, xyz[:, np.newaxis])
    return b[:, 0]


def grow_particle(pos1, pos2, bond, angle):
    '''
    Generate a random coordinates based on the coordinates of two parents, bond and angle

    https://scicomp.stackexchange.com/questions/27965/rotate-a-vector-by-a-randomly-oriented-angle
    @tparker

    Parameters
    ----------
    pos1 : np.ndarray
        coordinates of far parent
    pos2 : np.ndarray
        coordinates of near parent
    bond : float
        bond between near parent and new particle
    angle : float
        angle between far parent, near parent and new particle

    Returns
    -------
    xyz_new : array
        coordinates of new particle
    '''
    delta = pos2 - pos1
    v = delta / np.sqrt(np.dot(delta, delta))

    a = np.array([1, 0, 0])
    if np.dot(a, v) > 0.99:
        a = np.array([0, 1, 0])
    xb = np.cross(a, v)
    xb = xb / np.sqrt(np.dot(xb, xb))
    yb = np.cross(v, xb)

    phi = np.random.random() * 2 * np.pi
    w = np.sin(angle) * np.cos(phi) * xb + np.sin(angle) * np.sin(phi) * yb + np.cos(angle) * v

    return pos2 + w * bond


def periodic_distance(pos1, pos2, box, distance_max=None):
    '''
    Calculate the distance between two points under periodic boundary condition
    
    Parameters
    ----------
    pos1: np.ndarray
    pos2: np.ndarray
    box: np.ndarray of shape (3,)
        Lengths of rectangular periodic box. Triclinic box is not supported yet.
    distance_max : float, optional
        The maximum distance to be considered. Will return None if distance larger than it

    Returns
    -------
    distance : float or None
        May return None if distance is apparently larger than `distance_max`
    '''
    delta = pos1 - pos2
    if distance_max is not None and any((np.abs(delta) > distance_max) & (np.abs(delta) < box - distance_max)):
        return None

    ### elements of delta will be transformed to (-0.5, 0.5]
    delta -= np.ceil(delta / box - 0.5) * box

    ### elements of delta will be transformed to [-0.5, 0.5)
    # delta -= np.floor(delta / box + 0.5) * box

    if distance_max is not None and any(np.abs(delta) > distance_max):
        return None

    return math.sqrt(np.dot(delta, delta))


def find_clusters(elements, func):
    '''
    Group elements into clusters

    Parameters
    ----------
    elements : Iterable
    func : func
        This function has two arguments. Both should be element from `elements`.
        Return True if two elements should be grouped into the same cluster.

    Returns
    -------
    clusters : list of list of int
        The clusters formed by elements. Each cluster is represented by a list of the index of elements.
    '''
    n_element = len(elements)
    flag = [0] * n_element
    n_cluster = 0
    clusters = []
    for i in range(n_element):
        if flag[i]:
            continue
        path = {i}
        expand = set()
        n_cluster += 1
        flag[i] = n_cluster
        for j in range(i + 1, n_element):
            if flag[j]:
                continue
            if func(elements[i], elements[j]):
                path.add(j)
                expand.add(j)
                flag[j] = n_cluster
        while expand:
            m = expand.pop()
            for n in range(i + 1, n_element):
                if flag[n]:
                    continue
                if func(elements[m], elements[n]):
                    path.add(n)
                    expand.add(n)
                    flag[n] = n_cluster
        clusters.append(list(path))

    return clusters


def find_clusters_consecutive(elements, func):
    '''
    Group elements into clusters. If element i and j are in the same group, all elements between i and j will also be put in the same group.

    Parameters
    ----------
    elements : Iterable
    func : func
        This function has two arguments. Both should be element from `elements`.
        Return True if two elements should be grouped into the same cluster.

    Returns
    -------
    clusters : list of list of int
        The clusters formed by elements. Each cluster is represented by a list of the index of elements.
    '''
    n_element = len(elements)
    flag = [0] * n_element
    n_cluster = 0
    clusters = []
    for i in range(n_element):
        if not flag[i]:
            n_cluster += 1
            flag[i] = n_cluster
            clusters.append({i})
        for j in reversed(range(i + 1, n_element)):
            if flag[j]:
                break
            if func(elements[i], elements[j]):
                for n in range(i + 1, j + 1):
                    flag[n] = flag[i]
                    clusters[-1].add(n)
                break

    clusters = [list(c) for c in clusters]

    return clusters