"""
    ./py/dos.py

    Author: Fabian R. Lux
    Date:   2023-02-28

    Given a numpy array of eigenvalues, compute the density of states
"""

import numpy as np


def gaussian(x, s, mu):
    return np.exp(-0.5*(x-mu)**2/s**2) / np.sqrt(2*np.pi*s**2)

def density_of_states(spec, E_min, E_max, n_E=100, gamma=0.1):
    """
        spec: 1D numpy array containing the eigenvalues
    """

    # E_min = np.amin(spec)
    # E_max = np.amax(spec)

    E_mesh = np.linspace(E_min, E_max, n_E)


    dos = np.zeros(n_E)

    for i in range(n_E):
        mu = E_mesh[i]

        for E in spec:
            dos[i] += gaussian(E, gamma, mu)


    return E_mesh, dos