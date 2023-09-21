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

n = 1
k = 1

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

def spectrum(n,k):
    """
        exact diagonalization
    """

    mag = magnetic_representation(generators, k,n)

    H = scipy.sparse.csr_matrix((k*d, k*d), dtype=complex)

    for i in range(4*g):
        H += (-mag[i]) / (4*g)

    H_dense = H.todense()
    eigenvalues, eigenvectors = scipy.linalg.eigh(H_dense)

    return eigenvalues, eigenvectors


eigenvalues, eigenvectors = spectrum(n,k)

np.save("./out/tffg_3_localization_eigenvalues.npy",  eigenvalues)
np.save("./out/tffg_3_localization_eigenvectors.npy", eigenvectors)

