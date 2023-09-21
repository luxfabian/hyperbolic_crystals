"""
    ./py/post_process_tffg_localization.py

    Author: Fabian R. Lux
    Date:   9/11/2023
"""

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import RectBivariateSpline
from scipy.sparse import csr_matrix
import hyperbolic_disk as H

import tffg_iomodule as iomodule

from tffg_group_extension import get_extension,lil_matrix

from numba import njit
import scipy

n = 1
k = 7
g = 2

sigmas = get_extension(g,k,n)


access_point = iomodule.get_access_point()


d = access_point['d']

gammas = H.build_gammas(g)

fname_prefix = access_point['fprefix']

words_file = open(fname_prefix + ".words","r")

words = np.zeros((d, 2, 2),dtype=complex)

words_string = np.empty(d, dtype=str)
s_cocycle  = np.zeros((d,k,k),dtype=complex)

while True:

    
    # Get next line from file
    line = words_file.readline()

    # if line is empty
    # end of file is reached
    if not line:
        break

    rep = [int(n) for n in line.rstrip().split(' ')]

    word = np.eye(2,dtype=complex)

    word_string = ""

    s  = np.eye(k,dtype=complex)
    # print(line.split(' '))
    # reg_data = np.genfromtxt(line, dtype=int, delimiter=" ")

    for i in range(len(rep)-1):

        j = rep[i+1]

        if j>0:
            j = j-1
            word = word.dot(gammas[j])
            s = s.dot(sigmas[j].toarray())

            word_string += str(j)


    # print(rep)
    words[rep[0]] = word

    s_cocycle[rep[0]] = s

    words_string[rep[0]] = word_string
        
        
        # print("Line{}: {}".format(count, line.strip()))
    
words_file.close()


generators = []

for i in range(4*g):

    H = lil_matrix((d, d))

    reg_fname = fname_prefix + "_" + str(i) + ".reg"
    reg_data = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i, j] in reg_data:
        if i >= 0 and j >= 0:
            H[i, j] += 1


    generators.append(H)

generators =np.roll(generators, 2*g)

def isDiag(M):
    i, j = np.nonzero(M)
    return np.all(i == j)

magnetic_generators = []

for g in range(4*g):

    omega = np.zeros((d,d), dtype=complex)

    nz_i, nz_j = generators[g].nonzero()

    # for nz in range(d):

    #     h_m = nz_i[nz] # h_m = h_j g^-1 = \pi_R(g)_{mj} h_j

    #     i = nz_i[nz]
    #     j = nz_j[nz]

    #     #\omega(g_i g^-1, g)

    #     omega_ij = s_cocycle[i]
        

    for h_j in range(d):

        for jj in range(len(nz_j)):
            if h_j == nz_j[jj]:
                h_m = nz_i[jj] # h_m = h_j g^-1 = \pi_R(g)_{mj} h_j
                 
                omega_ij = s_cocycle[h_j]
                omega_ij = omega_ij.dot(np.linalg.inv(s_cocycle[g+1]))
                omega_ij = omega_ij.dot(np.linalg.inv(s_cocycle[h_m])) 
                
                if not isDiag(omega_ij):
                    #print("nondiagonal at", words_string[h_j])
                    print("break")
                    break
                # else:

                #     if np.abs(omega_ij[0,0] - 1) > 1e-5:

                #         print( np.diagonal(omega_ij))

                omega[:,h_j] = np.ones(d,dtype=complex) * omega_ij[0,0]
        

        # omega_ij = s_cocycle[m]
        # omega_ij = omega_ij.dot(s_cocycle[g+1])
        # omega_ij = omega_ij.dot(s_cocycle[j].transpose())

        



# H = scipy.sparse.csr_matrix((d, d), dtype=complex)


#     z0 = 0 + 0j

#     zs = np.zeros(d, dtype=complex)

#     for i in range(d):
#         zs[i] = H.moebius(words[i], z0)

#     eigenvalues = np.load("./out/tffg_"+str(mod)+"_localization_eigenvalues.npy")
#     eigenvectors = np.load("./out/tffg_"+str(mod)+"_localization_eigenvectors.npy")


#     def index_localization(i):

#         v = eigenvectors[:,i]
#         weights = v*v.conj()

#         z_mean = np.sum( zs * weights) / np.sum(weights)
#         z_var =  np.sqrt( np.sum( np.abs(zs - z_mean)**2 * weights) / np.sum(weights) )

#         return H.hyperbolic_distance(0,z_var)


#     localization = np.array( [index_localization(i) for i in range(d) ] )

#     plt.stairs(*np.histogram(localization, bins=bins, density=True), label=label)


# plot(3,20,"$p^n=3^1$")
# plot(5,200,"$p^n=5^1$")

# plt.legend()
# ax = plt.gca()
# # ax.set_xlabel(r"$d(0, \sum_i  z_i  \| \psi_i \|^2 ) $")

# ax.set_xlabel(r"$d(0,\Delta  z)$ with  $ \Delta z= \sqrt{ \langle | z - \langle z \rangle |^2 \rangle}$")
# ax.set_ylabel("Frequency")

# # plt.hist(localization)
# plt.grid()
# plt.tight_layout()
# plt.savefig('./out/tffg_localization.png',dpi=300)
# plt.clf() 