"""
    ./py/run_tffg_hofstaedter_kpm.py

    Author: Fabian R. Lux
    Date:   9/11/2023
"""
import tffg_iomodule as iomodule
import numpy as np

import scipy
from scipy.sparse import lil_matrix

from tffg_group_extension import magnetic_representation

import kpm

from math import floor

# -- parameters

k = 101
ns = np.array([n for n in range(floor(k/2))])

n_energies = 2048
n_moments = 512
n_random_states = 10

# -- group information

access_point = iomodule.get_access_point()

g = access_point['g']
N = access_point['N']
d = access_point['d']

fname_prefix = access_point['fprefix']


print("dimension", d)
print(fname_prefix)

# -- read generators

generators = []

couplings = np.array( [1,1,1,1], dtype=complex)
couplings = np.hstack( (couplings, couplings.conjugate()))

for i in range(4*g):

    H = lil_matrix((d, d))

    reg_fname = fname_prefix + "_" + str(i) + ".reg"
    reg_data = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i, j] in reg_data:
        if i >= 0 and j >= 0:
            H[i, j] += 1

    generators.append(H)


def spectrum(n):
    """
        exact diagonalization
    """

    mag = magnetic_representation(generators,k,n)

    H = scipy.sparse.csr_matrix((k*d, k*d), dtype=complex)

    # U = mag[0]
    # V = mag[1]

    # H = (-U.dot(U) - V.dot(V) + 16 * (U+V))/12 / 4

    # H += H.conjugate().transpose()

    for i in range(4*g):
        H += (-mag[i].dot(mag[i]) + 16 * mag[i])/12/(4*g)     

    emesh, dos = kpm.density_of_states(
            H, scale=2, n_moments=n_moments, n_energies=n_energies, n_random_states=n_random_states)

    return emesh, dos


hofstadter = np.zeros((2*len(ns), n_energies),dtype=float)
flux = np.zeros(2*len(ns),dtype=float)

for i in range(len(ns)):
    # -- flux
    phi = ns[i]/k

    print("ϕ={}/{}".format(ns[i],k))

    flux[i] = phi

    emesh, dos = spectrum(ns[i])
    hofstadter[i,:] = dos

    flux[2*len(ns)-i-1] = 1-phi
    hofstadter[2*len(ns)-i-1,:] = dos

# print(flux)

np.save("./out/hofstadter_kpm_emesh.npy", emesh)
np.save("./out/hofstadter_kpm_flux.npy", flux)
np.save("./out/hofstadter_kpm.npy", hofstadter)

