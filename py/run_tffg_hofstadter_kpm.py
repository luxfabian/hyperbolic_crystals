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

# -- parameters

k = 197
ns = np.array([n+1 for n in range(round((k+1)/2))])

n_energies = 500
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
for i in range(4*g):

    H = lil_matrix((d, d))

    reg_fname = fname_prefix + "_" + str(i) + ".reg"
    reg_data = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i, j] in reg_data:
        H[i, j] += 1

    generators.append(H)

def spectrum(n):
    """
        exact diagonalization
    """

    mag = magnetic_representation(generators, k,n)

    H = scipy.sparse.csr_matrix((k*d, k*d), dtype=complex)

    for i in range(4*g):
        H += (-mag[i])

    emesh, dos = kpm.density_of_states(
            H, scale=10, n_moments=n_moments, n_energies=n_energies, n_random_states=n_random_states)

    return emesh, dos


hofstadter = np.zeros((len(ns), n_energies),dtype=float)
flux = np.zeros(len(ns),dtype=float)

for i in range(len(ns)):
    # -- flux
    phi = ns[i]/k

    print("ϕ={}/{}".format(ns[i],k))

    flux[i] = phi

    emesh, dos = spectrum(ns[i])
    hofstadter[i,:] = dos


np.save("./out/hofstadter_kpm_emesh.npy", emesh)
np.save("./out/hofstadter_kpm_flux.npy", flux)
np.save("./out/hofstadter_kpm.npy", hofstadter)

