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
        position of far parent
    pos2 : np.ndarray
        position of near parent
    bond : float
        bond between near parent and new particle
    angle : float
        angle between far parent, near parent and new particle

    Returns
    -------
    pos3 : array
        position of new particle
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
    theta = np.pi - angle
    w = np.sin(theta) * np.cos(phi) * xb + np.sin(theta) * np.sin(phi) * yb + np.cos(theta) * v

    return pos2 + w * bond


def relocate_hydrogen(hydrogen):
    '''
    Update the position of a hydrogen atom according to the parent atom and other neighbors of parent

    Parameters
    ----------
    hydrogen : Atom
    '''
    parent = hydrogen.bond_partners[0]
    if len(parent.bonds) == 1:
        vector = np.random.random(3) - 0.5
        vector = vector / np.sqrt(np.dot(vector, vector))
        hydrogen.position = parent.position + vector * 0.1  # bXH = 0.1 nm
    elif len(parent.bonds) == 2:
        far = next(p for p in parent.bond_partners if p is not hydrogen)
        hydrogen.position = grow_particle(far.position, parent.position, 0.1,
                                          math.pi * 2 / 3)  # bXH = 0.1 nm, aYXH = 120 degree
    else:
        vector = np.array([0., 0., 0.])
        for atom in parent.bond_partners:
            if atom is not hydrogen:
                vector += atom.position - parent.position
        hydrogen.position = parent.position - vector / np.sqrt(np.dot(vector, vector)) * 0.1  # bXH = 0.1 nm
        hydrogen.position = hydrogen.position + (np.random.random(3) - 0.5) * 0.02  # add random noise up to 0.02 nm


def periodic_distance(pos1, pos2, box, distance_max=None):
    '''
    Calculate the distance between two points under periodic boundary condition
    
    Parameters
    ----------
    pos1 : np.ndarray
    pos2 : np.ndarray
    box : np.ndarray of shape (3,)
        Lengths of rectangular periodic box. Triclinic box is not supported yet.
    distance_max : float, optional
        The maximum distance to be considered. Will return None if distance larger than it

    Returns
    -------
    distance : float or None
        May return None if distance is apparently larger than `distance_max`
    '''
    delta = pos2 - pos1
    abs_delta = np.abs(delta)
    if distance_max is not None and any((abs_delta > distance_max) & (abs_delta < box - distance_max)):
        return None

    ### elements of delta will be transformed to (-0.5, 0.5]
    delta -= np.ceil(delta / box - 0.5) * box

    ### elements of delta will be transformed to [-0.5, 0.5)
    # delta -= np.floor(delta / box + 0.5) * box

    return math.sqrt(np.dot(delta, delta))


def periodic_distances(positions1, positions2, box):
    '''
    Calculate the distances between two groups of points under periodic boundary condition

    The length of positions1 and positions2 must be the same.
    The output will be the distances between point pairs with the same length as positions1 and positions2.
    This function makes use of the parallelization of numpy, and is preferred to calling `periodic_distance` in a loop.

    Parameters
    ----------
    positions1 : np.ndarray
    positions2 : np.ndarray
    box : np.ndarray of shape (3,)
        Lengths of rectangular periodic box. Triclinic box is not supported yet.

    Returns
    -------
    distances : np.ndarray
        The distances between each point pair
    '''
    delta_vectors = positions2 - positions1
    delta_vectors -= np.ceil(delta_vectors / box - 0.5) * box

    return np.sqrt(np.sum(delta_vectors ** 2, axis=1))


def periodic_angle(pos1, pos2, pos3, box):
    '''
    Calculate the angle between three points under periodic boundary condition

    Parameters
    ----------
    pos1 : np.ndarray
    pos2 : np.ndarray
    pos3 : np.ndarray
    box : np.ndarray of shape (3,)
        Lengths of rectangular periodic box. Triclinic box is not supported yet.

    Returns
    -------
    angle : float
    '''
    vec1 = pos1 - pos2
    vec2 = pos3 - pos2

    ### elements of delta will be transformed to (-0.5, 0.5]
    vec1 -= np.ceil(vec1 / box - 0.5) * box
    vec2 -= np.ceil(vec2 / box - 0.5) * box

    cos = vec1.dot(vec2) / np.sqrt(vec1.dot(vec1) * vec2.dot(vec2))
    return float(np.arccos(np.clip(cos, -1, 1)))


def periodic_dihedral(pos1, pos2, pos3, pos4, box):
    '''
    Calculate the dihedral between four points under periodic boundary condition

    Parameters
    ----------
    pos1 : np.ndarray
    pos2 : np.ndarray
    pos3 : np.ndarray
    pos4 : np.ndarray
    box : np.ndarray of shape (3,)
        Lengths of rectangular periodic box. Triclinic box is not supported yet.

    Returns
    -------
    angle : float
    '''
    vec1 = pos2 - pos1
    vec2 = pos3 - pos2
    vec3 = pos4 - pos3

    ### elements of delta will be transformed to (-0.5, 0.5]
    vec1 -= np.ceil(vec1 / box - 0.5) * box
    vec2 -= np.ceil(vec2 / box - 0.5) * box
    vec3 -= np.ceil(vec3 / box - 0.5) * box

    n1 = np.cross(vec1, vec2)
    n2 = np.cross(vec2, vec3)
    cos = n1.dot(n2) / np.sqrt(n1.dot(n1) * n2.dot(n2))
    value = float(np.arccos(np.clip(cos, -1, 1)))
    sign = 1 if vec1.dot(n2) >= 0 else -1
    return sign * value


def rotate_points(positions, reference, target):
    '''
    Rotate a group of points as a whole so that a reference vector aligns with a target vector

    https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d

    Parameters
    ----------
    positions : array
        The positions of the points
    reference : array
        The reference vector for the roratition
    target : array
        The target vector to align with

    Returns
    -------
    positions_rotated : array
        The positions after the rotation
    '''
    # unit vectors
    vec_unit = reference / np.linalg.norm(reference)
    tar_unit = target / np.linalg.norm(target)

    # determine the axis of rotation and angle of rotation
    axis = np.cross(vec_unit, tar_unit)
    sin = np.linalg.norm(axis)
    cos = np.dot(vec_unit, tar_unit)

    # construct the rotation matrix
    skew_matrix = np.array([[0, -axis[2], axis[1]],
                            [axis[2], 0, -axis[0]],
                            [-axis[1], axis[0], 0]])
    rotation_matrix = np.eye(3) + skew_matrix + (1 - cos) / sin ** 2 * np.dot(skew_matrix, skew_matrix)

    return np.array([np.dot(rotation_matrix, p) for p in positions])


def find_clusters(elements, func):
    '''
    Group elements into clusters

    This method is slow. Use :func:find_clusters_in_graph if the connectivity between elements is known.

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
    Group elements into clusters

    This method differs from :func:find_clusters in that it groups consecutive elements into the same cluster.
    If element i and j are in the same cluster, all elements between i and j will also be put in the same cluster.

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


def find_clusters_in_graph(n_node, edges):
    """
    Find clusters in a graph

    Parameters
    ----------
    n_node : int
        Number of nodes in the graph
    edges : list of set of int
        The edges of the graph. edges[i] is the set of nodes that are connected to node i.

    Returns
    -------
    clusters : list of list of int
        The clusters formed by nodes. Each cluster is represented by a list of the index of nodes.
    """
    flag = [False] * n_node
    clusters = []
    for i in range(n_node):
        if flag[i]:
            continue
        cluster = []
        stack = {i}
        while stack:
            node = stack.pop()
            if flag[node]:
                continue
            flag[node] = True
            cluster.append(node)
            stack.update(edges[node])
        clusters.append(cluster)

    return clusters
