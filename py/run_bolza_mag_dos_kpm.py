"""
    ./py/run_tffg_hofstaedter_kpm.py

    Author: Fabian R. Lux
    Date:   9/11/2023
"""
import tffg_iomodule as iomodule
import numpy as np

import scipy

from scipy.sparse import diags, lil_matrix, kron,eye

from tffg_group_extension import magnetic_representation

import kpm

import bolza_hoppings

# -- hopping matrices
TA0, TB0, TA1 = bolza_hoppings.bolza_hoppings()

# -- magnetic parameters
n = 2
k = 5

# -- KPM parameters
n_energies = 2048
n_moments = 512
n_random_states = 10

# -- group information

access_point = iomodule.get_access_point()

g = access_point['g']
N = access_point['N']
d = access_point['d']

if(g!=2):
    raise ValueError("This script only works for genus 2, however g={} is selected".format(g))

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

    H = scipy.sparse.csr_matrix((k*d*48, k*d*48), dtype=complex)
    e = eye(k*d)

    sA = 1 # np.exp(2*np.pi*1j * n/k/8)
    sB = 1 # np.exp(2*np.pi*1j * n/k/3)

    for i in range(8):
        H += - sA*kron(mag[i],TA1[i]) / 4

    H += - kron(e ,sA*TA0+sB*TB0) / 4

    H = H + H.conjugate().transpose()

    emesh, dos = kpm.density_of_states(
            H, scale=1.2, n_moments=n_moments, n_energies=n_energies, n_random_states=n_random_states)

    return emesh, dos


emesh, dos = spectrum(n)

np.save("./out/bolza_"+str(N)+"_"+str(n)+"_"+str(k)+"_dos_kpm_emesh.npy", emesh)
np.save("./out/bolza_"+str(N)+"_"+str(n)+"_"+str(k)+"_dos_kpm_flux.npy", n/k)
np.save("./out/bolza_"+str(N)+"_"+str(n)+"_"+str(k)+"_dos_kpm.npy", dos)

