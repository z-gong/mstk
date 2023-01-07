import math
import numpy as np
from mstk.topology import Atom

__all__ = [
    'calc_weighted_average',
    'calc_com',
    'calc_rg',
    'calc_hull_volume',
]


def calc_weighted_average(array, weight):
    '''
    Calculate the weighted average of 2-dimensional data

    Parameters
    ----------
    array : np.array
    weight : np.array
        If weight for all elements equal to zero. Will ignore weight

    Returns
    -------
    average : np.array
    '''
    if all(weight == 0):
        weight = np.ones(weight.shape)
    return np.sum(array * weight[:, np.newaxis], axis=0) / np.sum(weight)


def calc_com(atoms):
    '''
    Calculate the center of mass of a group of atoms

    Parameters
    ----------
    atoms : list of Atom

    Returns
    -------
    com : np.array
    '''
    positions = np.array([atom.position for atom in atoms])
    masses = np.array([atom.mass for atom in atoms])

    return calc_weighted_average(positions, masses)


def calc_rg(atoms):
    '''
    Calculate the radius of gyration of a group of atoms

    Parameters
    ----------
    atoms : list of Atom

    Returns
    -------
    rg : float
    '''
    positions = np.array([atom.position for atom in atoms])
    masses = np.array([atom.mass for atom in atoms])
    if all(masses == 0):
        masses.fill(1.0)

    com = calc_weighted_average(positions, masses)
    positions -= com
    rsq_ew = np.sum(positions ** 2, axis=1)

    return math.sqrt(sum(masses * rsq_ew) / sum(masses))


def calc_hull_volume(atoms):
    '''
    Calculate the volume occupied by a groups of at least four atoms

    Parameters
    ----------
    atoms : list of Atom

    Returns
    -------
    volume : float
    '''
    if len(atoms) < 4:
        raise Exception('At least 4 atoms required for 3D convex hull')

    from scipy.spatial import ConvexHull

    positions = np.array([atom.position for atom in atoms])
    hull = ConvexHull(positions)

    return float(hull.volume)
