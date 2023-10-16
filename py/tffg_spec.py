"""
    ./py/run_tffg_spec.py

    Author: Fabian R. Lux
    Date:   9/24/2023
"""
import tffg_iomodule as iomodule
import numpy as np

import scipy
from scipy.sparse import lil_matrix


from tffg_group_extension import magnetic_representation


def spec_mag_nn_ed(n=1, k=1):
    """
        Spectrum for magnetic representation and random next-nearest coupling
        obtained via exact diagonalization.

        flux = n/k
    """

    access_point = iomodule.get_access_point()

    g = access_point['g']
    N = access_point['N']
    d = access_point['d']

    fname_prefix = access_point['fprefix']

    # -- load generators from file
    generators = []
    for i in range(4*g):

        H = lil_matrix((d, d))

        reg_fname = fname_prefix + "_" + str(i) + ".reg"
        reg_data = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

        for [i, j] in reg_data:
            H[i, j] += 1

        generators.append(H)

    generators = np.roll(generators, 2*g)

    # -- construct magnetic representation

    mag = magnetic_representation(generators, k, n)

    H = scipy.sparse.csr_matrix((k*d, k*d), dtype=complex)

    couplings = []

    for i in range(4*g):
        H += -mag[i] / (4*g)

        # -- next-nearest neighbors
        for j in range(4*g):
            if abs(i-j) != 2*g:  # exclude inverse operation
                t_nn = np.random.rand()
                couplings.append(t_nn)

                H += -t_nn * mag[i].dot(mag[j]) / (4*g )

            for k in range(4*g):
                if abs(j-k) != 2*g:  # exclude inverse operation
                    t_nn = np.random.rand()
                    couplings.append(t_nn)

                    H += -t_nn * mag[i].dot(mag[j]).dot(mag[k]) / (4*g )

    H_dense = H.todense()
    eig = np.linalg.eigvalsh(H_dense)

    return eig

if __name__=='__main__':

    eig = spec_mag_nn_ed(1,3)

    np.save("./out/tffg_spec_ed.npy", eig)