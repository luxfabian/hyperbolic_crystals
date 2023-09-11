"""
    ./py/run_tffg_mag.py

    Author: Fabian R. Lux
    Date:   9/11/2023
"""
import tffg_iomodule as iomodule
import numpy as np

import scipy
from scipy.sparse import lil_matrix


from tffg_group_extension import magnetic_representation


# -- parameters
n = 1
k = 3


access_point = iomodule.get_access_point()

g = access_point['g']
N = access_point['N']
d = access_point['d']



fname_prefix = access_point['fprefix']

print("dimension", d)
print(fname_prefix)

generators = []
for i in range(4*g):

    H = lil_matrix((d, d))

    reg_fname = fname_prefix + "_" + str(i) + ".reg"
    reg_data = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i, j] in reg_data:
        H[i, j] += 1

    generators.append(H)


mag = magnetic_representation(generators, k,n)

H = scipy.sparse.csr_matrix((k*d, k*d), dtype=complex)

for i in range(4*g):
    H += mag[i]

H_dense = H.todense()
eig = scipy.linalg.eigvalsh(H_dense)

np.save("./out/tffg_spec_"+str(abs(N))+"_"+str(n)+"_"+str(k)+".npy", eig)

# # energy, dos = kpm.density_of_states(
# #     H, scale=h_scale, n_moments=n_moments, n_energies=n_energies, n_random_states=n_random_states
# # )

# # np.save("./out/hyperbolic_chern_energy.npy", energy)
# # np.save("./out/hyperbolic_chern_dos.npy", dos)
