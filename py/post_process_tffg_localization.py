"""
    ./py/post_process_tffg_localization.py

    Author: Fabian R. Lux
    Date:   9/11/2023
"""

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import RectBivariateSpline

import hyperbolic_disk as H

import tffg_iomodule as iomodule

from numba import njit

g = 2


def plot(mod, bins,label):

    access_point = iomodule.get_specific_access_point(g, mod)

    d = access_point['d']

    gammas = H.build_gammas(g)

    fname_prefix = access_point['fprefix']

    words_file = open(fname_prefix + ".words","r")

    words = np.zeros((d, 2, 2),dtype=complex)

    print(fname_prefix)

    while True:

    
        # Get next line from file
        line = words_file.readline()
    
        # if line is empty
        # end of file is reached
        if not line:
            break

        rep = [int(n) for n in line.rstrip().split(' ')]

        word = np.eye(2,dtype=complex)
        # print(line.split(' '))
        # reg_data = np.genfromtxt(line, dtype=int, delimiter=" ")

        for i in range(len(rep)-1):

            k = rep[i+1]

            if k>0:
                word = word.dot(gammas[k-1])

        words[rep[0]] = word
        
        
        # print("Line{}: {}".format(count, line.strip()))
    
    words_file.close()

    z0 = 0 + 0j

    zs = np.zeros(d, dtype=complex)

    for i in range(d):
        zs[i] = H.moebius(words[i], z0)

    eigenvalues = np.load("./out/tffg_"+str(mod)+"_localization_eigenvalues.npy")
    eigenvectors = np.load("./out/tffg_"+str(mod)+"_localization_eigenvectors.npy")


    def index_localization(i):

        v = eigenvectors[:,i]
        weights = v*v.conj()

        z_mean = np.sum( zs * weights) / np.sum(weights)
        z_var =  np.sqrt( np.sum( np.abs(zs - z_mean)**2 * weights) / np.sum(weights) )

        return H.hyperbolic_distance(0,z_var)


    localization = np.array( [index_localization(i) for i in range(d) ] )

    plt.stairs(*np.histogram(localization, bins=bins, density=True), label=label)


plot(3,20,"$p^n=3^1$")
plot(5,200,"$p^n=5^1$")

plt.legend()
ax = plt.gca()
# ax.set_xlabel(r"$d(0, \sum_i  z_i  \| \psi_i \|^2 ) $")

ax.set_xlabel(r"$d(0,\Delta  z)$ with  $ \Delta z= \sqrt{ \langle | z - \langle z \rangle |^2 \rangle}$")
ax.set_ylabel("Frequency")

# plt.hist(localization)
plt.grid()
plt.tight_layout()
plt.savefig('./out/tffg_localization.png',dpi=300)
plt.clf() 