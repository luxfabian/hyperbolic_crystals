"""
    ./scripts/Y_junction.py

    Author: Fabian R. Lux
    Date:   2023-03-08

    Sets up and diagonalizes the Y-junction Hamiltonian. Only works with open boundary conditions as of now.
"""
import re
import numpy as np

import kpm
import iomodule
import hyperbolic_disk

import scipy
from scipy.sparse import lil_matrix


import matplotlib.pyplot as plt
import matplotlib.colors as colors

# -- options -----------------------------

n_couplings = 200
n_energies = 1000
n_moments = n_energies
n_random_states = 10

access_point = iomodule.get_access_point()

fname_prefix = access_point['fprefix']

p = access_point['p']
q = access_point['q']

A = hyperbolic_disk.get_A(p,q)
B = hyperbolic_disk.get_B(p,q)

alpha = 2*np.pi/p

# -- orientation of Y junction
phi = alpha /2

r0 = hyperbolic_disk.get_r0(p,q)

seed = (r0-0.12) * np.exp(1j*alpha/2)

# -- length scale of domain wall
l = 1

# -- import basis

d = access_point['d']

basis = np.empty(d, dtype=object)
with open(fname_prefix+".words", "r") as basis_file:
    for line in basis_file:
        no = int(re.sub("[^0-9]", "", line).strip())
        word = re.sub("[^A-Z]", "", line).strip()

        basis[no] = word

# -- convert basis to coordinates
zs = np.array( [ hyperbolic_disk.parse_word(word, A, B, seed) for word in basis ] )

# -- read regular representation

# -- initialize the matrices
g  = {}
g['A']  = lil_matrix((d,d))
g['iA'] = lil_matrix((d,d))
g['B']  = lil_matrix((d,d))
g['iB'] = lil_matrix((d,d))
g['AB'] = lil_matrix((d,d))

generators = [ "A", "iA", "B", "iB" ,"AB"]

fname_prefix = access_point['fprefix']
for k in generators:
    reg_fname = fname_prefix + "_" + k + ".reg" 
    reg_data  = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i,j] in reg_data:
        g[k][i,j] += 1.0

# -- Laplace operator
Delta = (g['A'] + g['iA'] + g['B'] + g['iB'])/4

# -- A projectors
lam = np.exp(-1j * 2*np.pi/ 5)
PA1 = lam *  ( g['A'] )
PA2 = lam *  ( PA1 @ g['A'] )
PA3 = lam *  ( PA2 @ g['A'] )
PA4 = lam *  ( PA3 @ g['A'] )
PA0 = lam *  ( PA4 @ g['A'] )

PA_1 = ( PA0 + PA1 + PA2 + PA3 + PA4) / 5

# -- B projector
lam = np.exp(-1j * 2*np.pi/ 4)
PB1 = lam *  ( g['B'] )
PB2 = lam *  ( PB1 @ g['B'] )
PB3 = lam *  ( PB2 @ g['B'] )
PB0 = lam *  ( PB3 @ g['B'] )

PB_1 = ( PB0 + PB1 + PB2 + PB3 ) / 4

# -- C projector
lam = np.exp(-1j * 2*np.pi/ 2)
PC1 = lam *  ( g['AB'] )
PC0 = lam *  ( PC1 @ g['AB'] )

PC_1 = ( PC0 + PC1 ) / 2

#--build local hamiltonian

epsilon = 0.8
id = scipy.sparse.identity(d)
H_1 = epsilon *(id-2*PA_1) +  (1-epsilon)*Delta
H_2 = epsilon *(id-2*PB_1) +  (1-epsilon)*Delta
H_3 = epsilon *(id-2*PC_1) +  (1-epsilon)*Delta

# -- list of nonzero elements

def nonzero_to_set(lil_mat):
    i, j = lil_mat.nonzero()

    nz = len(i)

    nzs = { (i[k], j[k]) for k in range(nz)}
    
    return nzs

H_1_nz = nonzero_to_set(H_1)
H_2_nz = nonzero_to_set(H_2)
H_3_nz = nonzero_to_set(H_3)

H_nz = ( H_1_nz.union(H_2_nz) ).union(H_3_nz)


# -- build global hamiltonian
H = lil_matrix((d,d),dtype=complex)
for (i,j) in H_nz:

    z = hyperbolic_disk.midpoint(zs[i], zs[j])

    chis = hyperbolic_disk.chi(z,l,phi)

    H[i,j] = H_1[i,j] * chis[0] + H_2[i,j] * chis[1] + H_3[i,j] * chis[2]
    # print(i,j,z, chis)


H_dense  = H.todense()

eigenvalues, eigenvectors = np.linalg.eigh(H_dense)

np.save("junction_eigenvalues.npy",eigenvalues)
np.save("junction_eigenvectors.npy",eigenvectors)

np.save("junction_zs.npy",zs)





# energy, dos = kpm.density_of_states(H, scale=25, n_moments=n_moments, n_energies=n_energies, n_random_states=n_random_states)

# fig = plt.figure()
# ax=plt.gca()
# ax.set_xlim((-2.0,2.0))
# # ax.set_xticks([0,0.25, 0.50, 0.75, 1.00])   

# plt.plot(energy,dos)

# # plt.imshow(
# #     dosmat.transpose(), 
# #     extent=(0,1, -1.2, 1.2), 
# #     origin='lower', aspect='auto',cmap='Blues', interpolation='spline16', 
# #     norm=colors.SymLogNorm(
# #         linthresh=0.1, linscale=0.6, vmin=0, vmax=np.amax(dosmat), base=10
# #         )
# #     ) 

# # plt.colorbar()

# plt.tight_layout()
# plt.savefig('junction_dos.png',dpi=200)
# plt.show()
# plt.clf()