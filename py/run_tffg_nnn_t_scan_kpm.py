"""
    ./py/run_tffg_t_scan_kpm.py

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

n = 1
k = 3
# ns = np.array([n for n in range(floor(k/2))])

nt = 200

n_energies = 2048
n_moments = 2048
n_random_states = 10

# -- group information

access_point = iomodule.get_access_point()

g = access_point['g']
N = access_point['N']
d = access_point['d']

ts = np.linspace(-2,2,nt) 

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
        H[i, j] += 1

    generators.append(H)


def spectrum(t):
    """
        exact diagonalization
    """

    mag = magnetic_representation(generators, k,n)

    H = scipy.sparse.csr_matrix((k*d, k*d), dtype=complex)

    for i in range(4*g):
        H += (-mag[i])/(4*g)   

        for j in range(4*g):  
            if abs(i-j) != 2*g:
                H += t*(-mag[i].dot(mag[j]))/(4*g) 

    emesh, dos = kpm.density_of_states(
            H, scale=15, n_moments=n_moments, n_energies=n_energies, n_random_states=n_random_states)

    
    return emesh, dos


hofstadter = np.zeros((nt, n_energies),dtype=float)
# flux = np.zeros(2*len(ns),dtype=float)

for i in range(nt):
    
    emesh, dos = spectrum(ts[i])

    
    hofstadter[i,:] = dos


# print(flux)

np.save("./out/nnn_t_scan_kpm_emesh.npy", emesh)
np.save("./out/nnn_t_scan_kpm_ts.npy", ts)
np.save("./out/nnn_t_scan_kpm.npy", hofstadter)

print(np.amin(emesh),np.amax(emesh))