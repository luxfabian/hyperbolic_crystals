import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors

def gaussian(x, s, mu):
    return np.exp(-0.5*(x-mu)**2/s**2) / np.sqrt(2*np.pi*s**2)

eigenvalues = np.load("junction_eigenvalues.npy")
eigenvectors = np.load("junction_eigenvectors.npy")

# -- find smallest eigenvalue
s = 0.10
mu = 0
vec = np.zeros( len(eigenvalues ) , dtype=complex)

for i in range(len(eigenvalues)):
    vec += gaussian( eigenvalues[i], s, mu) * np.abs(eigenvectors[:,i])**2

print(eigenvalues[i])

np.savetxt("real_vec.dat", vec.real)
np.savetxt("imag_vec.dat", vec.imag)

plt.plot(np.abs(vec))
plt.show()

# print(np.amin(np.abs(eigenvalues)))
# plt.plot(eigenvalues)
# plt.show()


# n_phi = 301
# n_r = 100

# phis = np.linspace(0, 2*np.pi, n_phi)
# rs = np.linspace(0,0.99999,n_r)

# ax = plt.subplot(111, polar=True)

# ctf = ax.contourf(phis,rs,chis[region], levels=400, cmap='viridis')