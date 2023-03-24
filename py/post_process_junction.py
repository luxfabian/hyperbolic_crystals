import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors

import hyperbolic_disk

def gaussian(x, s, mu):
    return np.exp(-0.5*(x-mu)**2/s**2) / np.sqrt(2*np.pi*s**2)

eigenvalues = np.load("junction_eigenvalues.npy")
eigenvectors = np.load("junction_eigenvectors.npy")
zs = np.load("junction_zs.npy")

# -- find smallest eigenvalue
s = 0.2
mu = 0
vec = np.zeros( len(eigenvalues ) , dtype=float)

for i in range(len(eigenvalues)):
    vec += gaussian( eigenvalues[i], s, mu) * np.abs(eigenvectors[:,i])**2

#print(eigenvalues[i])

vmax = np.amax(vec)
np.savetxt("ldos_abs_mag.dat", vec/vmax)

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