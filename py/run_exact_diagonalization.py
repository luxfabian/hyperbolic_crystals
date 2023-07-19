"""
    ./py/run_exact_diagonalization.py

    Author: Fabian R. Lux
    Date:   2023-03-22

    Exact diagonalization of the Laplace operator with numpy.
"""

import iomodule
import numpy as np

from scipy.sparse import lil_matrix


access_point = iomodule.get_access_point()

d = access_point['d']
fname_prefix = access_point['fprefix']

# -- initialize the matrices
g = {}
g['A'] = lil_matrix((d, d))
g['iA'] = lil_matrix((d, d))
g['B'] = lil_matrix((d, d))
g['iB'] = lil_matrix((d, d))
g['AB'] = lil_matrix((d, d))

generators = ["A", "iA", "B", "iB", "AB"]

# -- read regular representation
fname_prefix = access_point['fprefix']
for k in generators:
    reg_fname = fname_prefix + "_" + k + ".reg"
    reg_data = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i, j] in reg_data:
        if i >= 0 and j >= 0:
            g[k][i, j] += 1.0

# -- Laplace operator
Delta = (g['A'] + g['iA'] + g['B'] + g['iB'])/4
Delta_dense = Delta.todense()

eig = np.linalg.eigvalsh(Delta_dense)

p = access_point['p']
q = access_point['q']
N = -access_point['N']
np.save("./out/{}_{}_open_{}".format(p, q, N), eig)
