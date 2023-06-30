"""
    ./py/open_boundary_conditions.py

    Author: Fabian R. Lux
    Date:   2023-03-22

    Sets up the Laplace operator for open bc.
"""

import iomodule
import numpy as np

from scipy.sparse import lil_matrix



access_point = iomodule.get_access_point()


d = access_point['d']
fname_prefix = access_point['fprefix']


# -- read regular representation

# -- initialize the matrices
g = {}
g['A'] = lil_matrix((d, d))
g['iA'] = lil_matrix((d, d))
g['B'] = lil_matrix((d, d))
g['iB'] = lil_matrix((d, d))
g['AB'] = lil_matrix((d, d))

generators = ["A", "iA", "B", "iB", "AB"]

fname_prefix = access_point['fprefix']
for k in generators:
    reg_fname = fname_prefix + "_" + k + ".reg"
    reg_data = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i, j] in reg_data:
        if i >= 0 and j >= 0:
            g[k][i, j] += 1.0

# -- Laplace operator
Delta = (g['A'] + g['iA'] + g['B'] + g['iB'])/4

Delta_dense  = Delta.todense()


eig = np.linalg.eigvalsh(Delta_dense)

# eig = scipy.linalg.eigvalsh(H_dense)

# H_dense = np.asfortranarray(H.toarray())
# print("H constructed")

# eig, _, _ = scipy.linalg.lapack.zheev(H_dense, compute_v=0)

p = access_point['p']
q = access_point['q']
N = -access_point['N']
np.save("./out/{}_{}_open_{}".format(p,q,N), eig)    