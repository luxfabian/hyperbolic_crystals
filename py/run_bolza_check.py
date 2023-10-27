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
n = 1
k = 3


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



mag = magnetic_representation(generators, k,n)

e = eye(k*d)


s = np.exp(n*2*np.pi * 1j / k)

xs = np.linspace(0,1, 100)


for x in xs:
    for y in xs:
        sA = np.exp(1j * 2*np.pi * x)
        sB = np.exp(1j * 2*np.pi * y)

        A = scipy.sparse.csr_matrix((k*d*48, k*d*48), dtype=complex)
        B = scipy.sparse.csr_matrix((k*d*48, k*d*48), dtype=complex)

        for i in range(8):
            A += sA*kron(mag[i],TA1[i]) 

        A += sB*kron(e,TA0) 
        # B =  kron(e,sB*TB0) 

        A_loop = eye(48*k*d)
        B_loop = eye(48*k*d)
        AB_loop = eye(48*k*d)
        for i in range(8):
            A_loop = A_loop.dot(A)
        # for i in range(3):
        #     B_loop = B_loop.dot(B)


        # AB_loop = A.dot(B).dot(A).dot(B)

        residual = A_loop - eye(48*k*d) / s
        print(abs(residual).max())

# residual = B_loop - eye(48*k*d)
# print(residual.min(),residual.max())

# residual = AB_loop - eye(48*k*d)
# print(residual.min(),residual.max())

# # print(A_loop*sA)

# import matplotlib.pyplot as plt
# plt.imshow(A_loop.toarray().real, interpolation="bilinear")
# plt.colorbar()
# plt.savefig("A_loop.png")

# # print(s)