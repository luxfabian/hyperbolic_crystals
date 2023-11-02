"""
    ./py/run_area_cocycle.py

    Author: Fabian R. Lux
    Date:   2023-03-08

    Sets up and diagonalizes the Y-junction Hamiltonian. 
    Only works with open boundary conditions as of now.
"""
import re
import numpy as np

import iomodule
import hyperbolic_disk

import scipy
from scipy.sparse import lil_matrix

# -- options -----------------------------

import sys

# print 'Number of arguments:', len(sys.argv), 'arguments.'
# print 'Argument List:', str(sys.argv)


theta =0.46
print(theta)

access_point = iomodule.get_access_point()

fname_prefix = access_point['fprefix']

p = access_point['p']
q = access_point['q']

A = hyperbolic_disk.get_A(p, q)
B = hyperbolic_disk.get_B(p, q)
iA = np.linalg.inv(A)
iB = np.linalg.inv(B)

alpha = 2*np.pi/p
beta  = 2*np.pi/q

area = p * ( np.pi - (alpha +  beta ))

# theta = theta*area



# -- orientation of Y junction
phi = alpha / 2

r0 = hyperbolic_disk.get_r0(p, q)

seed = (r0-0.12) * np.exp(1j*alpha/2)

# -- length scale of domain wall
l = 1

# -- import basis

d = access_point['d']

basis = np.empty(d, dtype=object)
with open(fname_prefix+".words", "r") as basis_file:
    for line in basis_file:
        no = int(re.sub("[^0-9]", "", line).strip())
        word = re.sub("[^A-Z]", "", line).strip()

        basis[no] = word

# -- convert basis to coordinates
zs = np.array([hyperbolic_disk.parse_word(word, A, B, seed) for word in basis])

#-- read regular representation

# -- initialize the matrices
g = {}
g['A'] = lil_matrix((d, d), dtype=complex)
g['iA'] = lil_matrix((d, d), dtype=complex)
g['B'] = lil_matrix((d, d), dtype=complex)
g['iB'] = lil_matrix((d, d), dtype=complex)
g['AB'] = lil_matrix((d, d), dtype=complex)

generators = ["A", "iA", "B", "iB", "AB"]

fname_prefix = access_point['fprefix']
for k in generators:
    reg_fname = fname_prefix + "_" + k + ".reg"
    reg_data = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i, j] in reg_data:
        if i >= 0 and j >= 0:
            z1 = seed
            z2 = hyperbolic_disk.parse_word(basis[j], A, B, seed)
            z3 = hyperbolic_disk.parse_word(basis[i], A, B, seed)
            
            g[k][i, j] = hyperbolic_disk.area_cocycle(z1, z2, z3, theta)

# -- Hamiltonian
H = (g['A'] + g['B']  + g['iA'] + g['iB'] )/4

# H = (H + H.conjugate().transpose())

residual = np.amax(np.abs((H - H.conjugate().transpose())) ) 
hermitian = residual < 1e-6

if not hermitian:
    print("Oh oh: the hamiltonian is not hermitian! The residual is {}".format(residual))
    exit


H_dense = H.todense()

eigenvalues, eigenvectors = np.linalg.eigh(H_dense)

np.save("./out/mag_eigenvalues.npy", eigenvalues)
np.save("./out/mag_eigenvectors.npy", eigenvectors)
np.save("./out/mag_zs.npy", zs)
