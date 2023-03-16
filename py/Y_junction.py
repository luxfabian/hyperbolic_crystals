"""
    ./scripts/Y_junction.py

    Author: Fabian R. Lux
    Date:   2023-03-08

    Sets up and diagonalizes the Y-junction Hamiltonian. Only works with open boundary conditions as of now.
"""
import numpy as np

import hyperbolic_disk
import re

p = 5
q = 4

A = hyperbolic_disk.get_A(p,q)
B = hyperbolic_disk.get_B(p,q)

fprefix = "../data/5_4_open_3"

alpha = 2*np.pi/p

# -- orientation of Y junction
phi = alpha /2

seed = 0.1 * np.exp(1j*alpha/2)

# -- length scale of domain wall
l = 0.25

# -- import basis

info_file = open(fprefix+".info")
info_lines = info_file.readlines()
dim = int(re.sub("[^0-9]","", info_lines[0]).strip() )
print(dim)
info_file.close()

basis = np.empty(dim, dtype=object)
with open(fprefix+".words", "r") as basis_file:
    for line in basis_file:
        no = int(re.sub("[^0-9]", "", line).strip())
        word = re.sub("[^A-Z]", "", line).strip()

        basis[no] = word

# -- convert basis to coordinates
z = np.array( [ hyperbolic_disk.parse_word(word, A, B, seed) for word in basis ] )

# -- read regular representation
