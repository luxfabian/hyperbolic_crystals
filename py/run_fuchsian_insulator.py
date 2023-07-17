"""
    ./py/run_fuchsian_insulator.py

    Author: Fabian R. Lux
    Date:   6/28/2023

    We implement the model proposed in 

    https://doi.org/10.1038/s41467-023-36767-8

    for a hyperbolic Chern insulator with non-trivial second Chern number. We investigate,
    whether the band gaps of the system remain open in the thermodynamic limit.

"""
import fuchsian_iomodule as iomodule
import numpy as np

import scipy
from scipy.sparse import lil_matrix

import kpm

from clifford_algebra import gamma

# -- parameters

h_scale = 6
m = 0.7
a = 0.2

n_energies = 500
n_moments = 500
n_random_states = 20

access_point = iomodule.get_access_point()

d = access_point['d']
p = access_point['p']
n = access_point['n']
fname_prefix = access_point['fprefix']

generators = [1, 2, 3, 4, 5, 6, 7, 8]

S = {}
for g in generators:

    H = lil_matrix((d, d))

    reg_fname = fname_prefix + "_" + str(g) + ".reg"
    reg_data = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i, j] in reg_data:
        H[i, j] += 1

    S[g-1] = H


gamma14 = gamma[0].dot(gamma[3])

H = scipy.sparse.lil_matrix((4*d, 4*d), dtype=complex)
for i in range(4):
    H += scipy.sparse.kron(scipy.sparse.lil_matrix(gamma[i]), S[i] - S[i+4]) / (2*1j) \
        + scipy.sparse.kron(scipy.sparse.lil_matrix(gamma[4]),  S[i] + S[i+4]) / 2

H += m * scipy.sparse.lil_matrix(np.kron(gamma[4], np.eye(d, dtype=complex)))

H += 1j*a * scipy.sparse.lil_matrix(np.kron(gamma14, np.eye(d, dtype=complex)))

H_dense = H.todense()


eig = np.linalg.eigvalsh(H_dense)

# eig = scipy.linalg.eigvalsh(H_dense)

# H_dense = np.asfortranarray(H.toarray())
# print("H constructed")

# eig, _, _ = scipy.linalg.lapack.zheev(H_dense, compute_v=0)

np.save("./out/hyperbolic_chern_eig_"+str(p)+"_"+str(n)+".npy", eig)

# energy, dos = kpm.density_of_states(
#     H, scale=h_scale, n_moments=n_moments, n_energies=n_energies, n_random_states=n_random_states
# )

# np.save("./out/hyperbolic_chern_energy.npy", energy)
# np.save("./out/hyperbolic_chern_dos.npy", dos)
