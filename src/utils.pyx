import numpy as np
cimport numpy as np

from iteration_utilities import duplicates

cpdef np.ndarray get_mvs_from_many(np.ndarray mvs, int lin, int lout, int l):
    """
    Get specific multiple vectors given a list of many (output from mvs).

    Args:
        mvs (np.ndarray): array containing multipole vectors [theta, phi].
        lin (int): initial multipole.
        lout (int): final multipole.
        l (int): multipole.

    Returns:
        np.ndarray.
    """
    cdef np.ndarray multipoles = np.arange(lin, lout + 1, 1, dtype=np.int32)
    cdef np.ndarray indexes = np.repeat(multipoles, 2*multipoles)
    return mvs[np.where(indexes == l)]

cpdef np.ndarray to_cart(np.ndarray theta_phi_array):
    """
    Convert spherical to Cartesian coordinate. This function works only for a single multipole or stacked arrays, e.g. 2 and 3 MVs together.

    Args:
        theta_phi_array (np.ndarray): array containing in radians [theta, phi].

    Returns:
        np.ndarray.
    """
    cdef np.ndarray x = np.sin(theta_phi_array[:, 0]) * np.cos(theta_phi_array[:, 1])
    cdef np.ndarray y = np.sin(theta_phi_array[:, 0]) * np.sin(theta_phi_array[:, 1])
    cdef np.ndarray z = np.cos(theta_phi_array[:, 0])
    return np.vstack((x, y, z)).T

cpdef np.ndarray mvs_north(np.ndarray mvs):
    """
    Get multipole vectors on north hemisphere.

    Args:
        mvs (np.ndarray): array containing multipole vectors [theta, phi].

    Returns:
        np.ndarray [theta, phi] in radians.
    """
    cdef int half_length = len(mvs) // 2
    cdef np.ndarray mvs_n = mvs[mvs[:, 0] < np.pi / 2]

    if len(mvs_n) == half_length:
        return mvs_n

    if len(mvs_n) == half_length + 1:
        duplicate_thetas = duplicates(mvs_n[:, 0])
        if len(duplicate_thetas) == 1:
            idx_to_remove = np.random.choice(np.where(mvs_n[:, 0] == duplicate_thetas)[0], size=1)[0]
            mvs_n = np.delete(mvs_n, idx_to_remove, axis=0)
        else:
            idx_to_remove = np.argmax(mvs_n[:, 0])
            mvs_n = np.delete(mvs_n, idx_to_remove, axis=0)
        return mvs_n

    if len(mvs_n) == half_length - 1:
        mvs_s = mvs[mvs[:, 0] >= np.pi / 2]
        sorted_mvs_s = mvs_s[np.argsort(mvs_s[:, 0])][:2]
        duplicate_thetas = duplicates(sorted_mvs_s[:, 0])
        if len(duplicate_thetas) == 1:
            idx_to_add = np.random.choice([0, 1])
            mvs_n = np.vstack((mvs_n, sorted_mvs_s[idx_to_add]))
        else:
            mvs_n = np.vstack((mvs_n, sorted_mvs_s[0]))
        return mvs_n

    raise ValueError(f"Found a bug! len(mvs_n)={len(mvs_n)}, expected={half_length} or {half_length+1} or {half_length-1}")
