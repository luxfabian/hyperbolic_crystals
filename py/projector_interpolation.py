"""
    ./scripts/projector_interpolation.py

    Author: Fabian R. Lux
    Date:   2023-02-28

    Interpolates between to projection operators and calculates the spectrum along the path
"""

import iomodule
import kpm
import scipy
import numpy as np
from scipy.sparse import lil_matrix

import matplotlib.pyplot as plt


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
    reg_data  = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i,j] in reg_data:
        g[k][i,j] += 1.0


# -- Laplace operator
Delta = g['A'] + g['iA'] + g['B'] + g['iB'] + 2 * g['AB']


energy, dos = kpm.density_of_states(Delta, scale=7, n_moments=100, n_energies=1000, n_random_states=30)
plt.plot(energy, dos)
plt.show()

# spec = scipy.linalg.eigh(Delta.todense(), eigvals_only=True)

# plt.plot(spec)
# plt.show()
# H_fname = fname_prefix + ".hamiltonian"

# with open(H_fname, 'w') as H_file:
#     for i in range(d):
#         for j in H.rows[i]:
#             line = fortran_format(i) + " " + fortran_format(j) + " " + fortran_format(H[i,j]) + "\n"
#             H_file.write(line)