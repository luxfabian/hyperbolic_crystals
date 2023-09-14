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

# -- parameters

k = 51
ns = np.array([n for n in range(k+1)])

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

    H_dense = H.todense()
    spec = scipy.linalg.eigh(H_dense, eigvals_only=True)

    return spec


hofstadter = np.zeros((len(ns), k*d),dtype=float)
flux = np.zeros(len(ns),dtype=float)

for i in range(len(ns)):
    # -- flux
    phi = ns[i]/k

    print("ϕ={}/{}".format(ns[i],k))

    flux[i] = phi
    hofstadter[i,:] = spectrum(ns[i])


np.save("./out/hofstadter_flux.npy", flux)
np.save("./out/hofstadter.npy", hofstadter)

