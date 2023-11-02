"""
    ./py/dos.py

    Author: Fabian R. Lux
    Date:   2023-02-28

    Given a numpy array of eigenvalues, compute the density of states as sum of Gaussians.
"""

import numpy as np

from numba import njit

@njit
def gaussian(x, s, mu):
    return np.exp(-0.5*(x-mu)**2/s**2) / np.sqrt(2*np.pi*s**2)


def density_of_states(spec, E_min, E_max, n_E=100, gamma=0.1):
    """
        spec: 1D numpy array containing the eigenvalues
    """

    E_mesh = np.linspace(E_min, E_max, n_E)

    dos = np.zeros(n_E)

    for i in range(n_E):
        mu = E_mesh[i]

        for E in spec:
            dos[i] += gaussian(E, gamma, mu) / len(spec)

    return E_mesh, dos

# @njit
# def weighted_density_of_states(spec, E_min, E_max, weights, n_E=100, gamma=0.1):
#     """
#         spec: 1D numpy array containing the eigenvalues
#     """

#     E_mesh = np.linspace(E_min, E_max, n_E)

#     dos = np.zeros(n_E)

#     for i in range(n_E):
#         mu = E_mesh[i]

#         for j in range(spec):
#             dos[i] += weights[j]*gaussian(spec[j], gamma, mu) / len(spec)

#     return E_mesh, dos


@njit
def weighted_density_of_states(eigenvalues, eigenvectors, E_min, E_max, weights, n_E=100, gamma=0.1):
    """
        spec: 1D numpy array containing the eigenvalues
    """

    E_mesh = np.linspace(E_min, E_max, n_E)

    dos = np.zeros(n_E)

    for i in range(n_E):
        mu = E_mesh[i]

        for j in range(len(eigenvalues)):
            for k in range(len(eigenvectors)):
                dos[i] += weights[k]*gaussian( eigenvalues[j], gamma, mu) * np.abs(eigenvectors[k,j])**2 / len(eigenvalues)

    return E_mesh, dos


@njit
def local_density_of_states(eigenvalues, eigenvectors, site, E_min, E_max, n_E=100, gamma=0.1):
    """
        spec: 1D numpy array containing the eigenvalues
    """

    E_mesh = np.linspace(E_min, E_max, n_E)

    dos = np.zeros(n_E)

    for i in range(n_E):
        mu = E_mesh[i]

        for j in range(len(eigenvalues)):
            dos[i] += gaussian( eigenvalues[j], gamma, mu) * np.abs(eigenvectors[site,j])**2 / len(eigenvalues)

    return E_mesh, dos