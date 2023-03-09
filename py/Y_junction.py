"""
    ./scripts/Y_junction.py

    Author: Fabian R. Lux
    Date:   2023-03-08

    Sets up - and diagonalizes the Y-junction Hamiltonian. Only works with open boundary conditions as of now.
"""
import numpy as np

import hyperbolic_disk

p = 5
q = 4

alpha = 2*np.pi/p

# -- orientation of Y junction
phi = alpha /2


# -- length scale of domain wall
l = 0.25

# -- import