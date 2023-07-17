"""
    ./py/projector_interpolation.py

    Author: Fabian R. Lux
    Date:   2023-02-28

    Interpolates between to projection operators and calculates the spectrum along the path
"""

# -- options -----------------------------

n_couplings = 200
n_energies = 200
n_moments = n_energies
n_random_states = 10

# ----------------------------------------

import iomodule
import kpm
import scipy
import numpy as np
from scipy.sparse import lil_matrix

import matplotlib.pyplot as plt
import matplotlib.colors as colors

plt.rcParams.update({'font.size': 20})

access_point = iomodule.get_access_point()

print("Projector interpolation")
print(access_point)

projectors = ["A", "B"]
orders     = [access_point['p'], access_point['q']]

generators = [ "A", "iA", "B", "iB" ,"AB"]

d = access_point['d']

# -- initialize the matrices
g  = {}
g['A']  = lil_matrix((d,d))
g['iA'] = lil_matrix((d,d))
g['B']  = lil_matrix((d,d))
g['iB'] = lil_matrix((d,d))
g['AB'] = lil_matrix((d,d))

fname_prefix = access_point['fprefix']
for k in generators:
    reg_fname = fname_prefix + "_" + k + ".reg"
    reg_data = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i, j] in reg_data:
        if i >= 0 and j >= 0:
            g[k][i, j] += 1.0



# -- Laplace operator
Delta = (g['A'] + g['iA'] + g['B'] + g['iB'])/4

id = scipy.sparse.identity(d)

# -- A projectors
lam = np.exp(-1j * 2*np.pi/ 5)
PA1 = lam *  ( g['A'] )
PA2 = lam *  ( PA1 @ g['A'] )
PA3 = lam *  ( PA2 @ g['A'] )
PA4 = lam *  ( PA3 @ g['A'] )
PA0 = lam *  ( PA4 @ g['A'] )

P = ( PA0 + PA1 + PA2 + PA3 + PA4) / 5 + 0.1 * Delta

P_dense  = P.todense() 

eig = np.linalg.eigvalsh(P_dense)

# eig = scipy.linalg.eigvalsh(H_dense)

# H_dense = np.asfortranarray(H.toarray())
# print("H constructed")

# eig, _, _ = scipy.linalg.lapack.zheev(H_dense, compute_v=0)

p = access_point['p']
q = access_point['q']
N = abs(access_point['N'])
np.save("./out/projector_{}_{}_open_{}".format(p,q,N), eig)    