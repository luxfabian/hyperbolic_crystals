"""
    ./py/run_tffg_hofstaedter_kpm.py

    Author: Fabian R. Lux
    Date:   9/11/2023
"""
import tffg_iomodule as iomodule
import numpy as np

import scipy
from scipy.sparse import lil_matrix, kron, eye

from tffg_group_extension import magnetic_representation

import kpm

from clifford_algebra import gamma

# -- parameters

t1 = 0.5
t2 = 0.2
n = 1
k = 3
# 2003/
# 33334/100003
n_energies = 2048
n_moments = 1024
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

    H = scipy.sparse.csr_matrix((k*d*4, k*d*4), dtype=complex)

    for i in range(2*g):
        H += kron(gamma[i], mag[i])/ (4*g) / 1j + kron(gamma[4], mag[i])/ (4*g) 

    H += H.conjugate().transpose()

    H +=  1 * kron(gamma[4],eye(k*d))
        # for j in range(4*g):  
        #     if abs(i-j) != 2*g:
        #         H += t1*(-mag[i].dot(mag[j]))/(4*g) 

        #         for m in range(4*g):  
        #             if abs(m-j) != 2*g:
        #                 H += t2*(-mag[i].dot(mag[j].dot(mag[m])))/(4*g) 

    emesh, dos = kpm.density_of_states(
            H, scale=3, n_moments=n_moments, n_energies=n_energies, n_random_states=n_random_states)

    return emesh, dos


emesh, dos = spectrum(n)

np.save("./out/tffg_nnn_"+str(g)+"_"+str(N)+"_"+str(k)+"_dos_kpm_emesh.npy", emesh)
np.save("./out/tffg_nnn_"+str(g)+"_"+str(N)+"_"+str(k)+"_dos_kpm_flux.npy", n/k)
np.save("./out/tffg_nnn_"+str(g)+"_"+str(N)+"_"+str(k)+"_dos_kpm.npy", dos)

