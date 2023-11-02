"""
    ./py/tffg_.py

    Author: Fabian R. Lux
    Date:   9/11/2023
"""
import tffg_iomodule as iomodule
import numpy as np

import scipy

from tffg_group_extension import magnetic_representation


import matplotlib.pyplot as plt


def operators(n=1,k=3):

    # -- group information

    access_point = iomodule.get_access_point()

    g = access_point['g']
    N = access_point['N']
    d = access_point['d']


    fname_prefix = access_point['fprefix']

    # -- read generators

    generators = []

    for i in range(4*g):

        H = scipy.sparse.lil_matrix((d, d))

        reg_fname = fname_prefix + "_" + str(i) + ".reg"
        reg_data = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

        for [i, j] in reg_data:
            H[i, j] += 1

        generators.append(H)

    generators = np.roll(generators, 2*g)

    mag = magnetic_representation(generators, k,n)

    # -- nearest neighbors
    H_1 = []
    H_1_mag = []
    
    # -- next-nearest neighbors
    H_2 = []
    H_2_mag = []

    for i in range(4*g): 
        H_1.append(  (-generators[i])/(4*g)   )
        H_1_mag.append(  (-mag[i])/(4*g)   )

        for j in range(4*g):  
            if abs(i-j) != 2*g:
                H_2.append( (-generators[i].dot(generators[j]))/(4*g) )
                H_2_mag.append( (-mag[i].dot(mag[j]))/(4*g) )

    return H_1, H_2, H_1_mag, H_2_mag


def spectrum(H_1, H_2, H_1_mag, H_2_mag):

    H     = scipy.sparse.csr_matrix(H_1[0], dtype=complex)
    H_mag = scipy.sparse.csr_matrix(H_1_mag[0], dtype=complex)

    t_1 = []
    t_2 = []

    for i in range(len(H_1)):

        t = np.random.rand() + 1j * np.random.rand()
        t_1.append(t)

        H += t * H_1[i]
        H_mag += t * H_1_mag[i]

    for i in range(len(H_2)):

        t = np.random.rand() + 1j * np.random.rand()
        t = 0.001*t
        t_2.append(t)

        H += t * H_2[i]
        H_mag += t * H_2_mag[i]

    # -- hermitianize
    H = (H + H.conj().transpose() ) / 2.0
    H_mag = (H_mag + H_mag.conj().transpose() ) / 2.0

    # -- diagonalize
    H_dense = H.todense()
    H_mag_dense = H_mag.todense()

    eig = scipy.linalg.eigvalsh(H_dense)
    eig_mag = scipy.linalg.eigvalsh(H_mag_dense)
    
    return np.sort(eig), np.sort(eig_mag)

def get_gaps(spec):

    bandwidth = np.amax(spec) - np.amin(spec)

    gaps = ( (np.roll(spec,-1) - spec )[0:-1:1] ) / bandwidth

    return gaps

if __name__=="__main__":

    H_1, H_2, H_1_mag, H_2_mag = operators()

    spec, spec_mag = spectrum(H_1, H_2, H_1_mag, H_2_mag)

    gaps = get_gaps(spec)
    gaps_mag = get_gaps(spec_mag)

    hist = np.histogram(gaps, bins=100, density=False)
    hist_mag = np.histogram(gaps_mag, bins=100, density=False)

    # plt.plot(np.sort(gaps), 'o')
   
    plt.stairs(*hist, label="$\phi=0$")
    plt.stairs(*hist_mag, label="$\phi=1/3$")

    ax = plt.gca()
    ax.set_axisbelow(True)

    plt.legend()
    plt.grid()

    ax.set_yscale("log")
    # ax.set_xlabel(r"$d(0, \sum_i  z_i  \| \psi_i \|^2 ) $")

    ax.set_xlabel("$\Delta E / (E_{max} - E_{min})$")
    ax.set_ylabel("Frequency")

    plt.tight_layout()
    plt.savefig("./out/tffg_gap_search_hist.png",dpi=300)
    plt.clf()

    