"""
    ./py/post_process_junction.py

    Author: Fabian R. Lux
    Date:   2023-03-08

    Analyzes the spectrum of the Y-junction calculation
"""
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors

from dos import density_of_states, weighted_density_of_states, local_density_of_states

import hyperbolic_disk

def gaussian(x, s, mu):
    return np.exp(-0.5*(x-mu)**2/s**2) / np.sqrt(2*np.pi*s**2)

eigenvalues = np.load("./out/mag_eigenvalues.npy")
eigenvectors = np.load("./out/mag_eigenvectors.npy")
zs = np.load("./out/mag_zs.npy")

ds = hyperbolic_disk.hyperbolic_distance(0, zs) 

# print(ds)

ell = 1

print(ds[0])

print(np.amin(ds))
weights = np.zeros(len(ds))

def sigmoid(d, ell):

    if d > ell:
        return 0
    else:
        return 1

    # return 1.0 / ( 1.0 + np.exp(+d-ell))

for i in range(len(ds)):

    weights[i] = sigmoid(ds[i], ell)

# print(weights)
# -- find smallest eigenvalue
s = 0.10
mu = 0
vec = np.zeros( len(eigenvalues ) , dtype=float)

# for i in range(len(eigenvalues)):
#     vec += gaussian( eigenvalues[i], s, mu) * np.abs(eigenvectors[:,i])**2

#print(eigenvalues[i])

vmax = np.amax(vec)
# np.savetxt("./out/ldos_abs_mag.dat", vec/vmax)


# -- plot the spectrum

# plt.plot(eigenvalues)
# ax = plt.gca()
# ax.set_ylim((-0.2,0.2))
# ax.set_xlim((3000,3250))
# plt.show()

E_mesh, dos    = density_of_states(eigenvalues, -1.2,1.2, n_E=500, gamma=0.005)

E_mesh, weighted_dos    = local_density_of_states(eigenvalues, eigenvectors, 0,-1.2,1.2, n_E=500, gamma=0.005)

# E_mesh, dos_10 = density_of_states(eigenvalues, -1.2,1.2, n_E=200, gamma=0.02)
# E_mesh, dos_20 = density_of_states(eigenvalues, -1.2,1.2, n_E=200, gamma=0.03)
# E_mesh, dos_30 = density_of_states(eigenvalues, -1.2,1.2, n_E=200, gamma=0.04)


plt.plot(E_mesh, dos/np.amax(dos), label=r'$s=0.01$')
plt.plot(E_mesh, weighted_dos/np.amax(weighted_dos), label=r'$s=0.01$ $\ell=0.1$')
# plt.plot(E_mesh, dos_10/np.amax(dos_10), label=r'$s=0.02$')
# plt.plot(E_mesh, dos_20/np.amax(dos_20), label=r'$s=0.03$')
# plt.plot(E_mesh, dos_30/np.amax(dos_30), label=r'$s=0.04$')
ax = plt.gca()
ax.set_xlabel(r"Energy")
ax.set_ylabel(r"Density of states")
# ax.set_ylim((0,1.2))
# ax.set_aspect(2)
#plt.plot(E_mesh, gauss_dat)
plt.legend(fontsize=10, loc='upper right')
plt.savefig("./out/mag_dos.png", dpi=300)
# plt.show()

# plt.plot(np.abs(vec))
# # plt.show()
# plt.clf()

# phis = np.zeros(len(eigenvalues))
# for i in range(len(eigenvalues)):
#     ldos =  vec[i]
#     z = zs[i]
#     r = abs(z)
#     phis[i] = hyperbolic_disk.arg(z)

# n_phi = 21
# phi_mesh = np.linspace(-np.pi, np.pi, n_phi, endpoint=False)

# dphi = phi_mesh[1]-phi_mesh[0]

# phi_mesh = phi_mesh + 0.02#+ dphi/2.0

# idos = np.zeros(n_phi)


# for i in range(n_phi):
#     idos[i] = 0.0

#     for j in range(len(eigenvalues)):
#         if phis[j] >= phi_mesh[i] and phis[j] < phi_mesh[i] + dphi:
#             idos[i] += vec[j]/vmax

# plt.plot(phi_mesh, idos)
# #plt.show()
# plt.clf()


# plt.bar(phi_mesh + dphi/2, idos, width=dphi)
# plt.show()

# # Compute pie slices
# N = 20
# theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
# radii = 10 * np.random.rand(N)
# # width = np.pi / 4 * np.random.rand(N)
# colors = plt.cm.Blues(idos)

# ax = plt.subplot(projection='polar')
# ax.bar(phi_mesh + dphi/2, idos, width=dphi, bottom=0.0, color=colors, alpha=0.5)

# deg = 180 / np.pi
# phi = np.pi/5
# p1, p2, p3  = hyperbolic_disk.junction_angles(phi)
# ax.set_thetagrids([p1*deg, p2*deg, p3*deg])
# ax.set_xticklabels([r"$\Upsilon_1^\infty$", r"$\Upsilon_2^\infty$", r"$\Upsilon_3^\infty$"], fontsize=20)
# ax.set_rticks([])
# ax.set_yticklabels([])


# # plt.colorbar(label="DOS")
# plt.tight_layout()
# plt.savefig("Y_junction_polar_dos.png", dpi=300)
# plt.show()

# -- local spectra

#idea: weight

# d_max = hyperbolic_disk.hyperbolic_distance(0,np.amax(np.abs(zs) ) )

# d_select = 3

# d_spread = 1



# print("Largest available distance:", d_max )

# print("Selected distance:", d_max )

# print("Spread: ", d_spread)

# n_pts = 10

# z_abs = #

# trajectors = d_s