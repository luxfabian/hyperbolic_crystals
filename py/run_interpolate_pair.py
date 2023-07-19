"""
    ./py/run_interpolate_pair.py

    Author: Fabian R. Lux
    Date:   2023-02-28

    Interpolates between two operators and calculates the spectrum along the path.
"""

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
orders = [access_point['p'], access_point['q']]

generators = ["A", "iA", "B", "iB", "AB"]

d = access_point['d']

# -- initialize the matrices
g = {}
g['A'] = lil_matrix((d, d))
g['iA'] = lil_matrix((d, d))
g['B'] = lil_matrix((d, d))
g['iB'] = lil_matrix((d, d))
g['AB'] = lil_matrix((d, d))

fname_prefix = access_point['fprefix']
for k in generators:
    reg_fname = fname_prefix + "_" + k + ".reg"
    reg_data = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i, j] in reg_data:
        g[k][i, j] += 1.0


# -- Laplace operator
Delta = g['A'] + g['iA'] + g['B'] + g['iB']  # + 2 * g['AB']

standard_coupling = 0.1
id = scipy.sparse.identity(d)

# -- A projector
lam = np.exp(-1j * 2*np.pi / 5)
PA1 = lam * (g['A'])
PA2 = lam * (PA1 @ g['A'])
PA3 = lam * (PA2 @ g['A'])
PA4 = lam * (PA3 @ g['A'])
PA0 = lam * (PA4 @ g['A'])

PA_1 = (PA0 + PA1 + PA2 + PA3 + PA4) / 5

lam = np.exp(-1j * 2 * 2*np.pi / 5)
PA1 = lam * (g['A'])
PA2 = lam * (PA1 @ g['A'])
PA3 = lam * (PA2 @ g['A'])
PA4 = lam * (PA3 @ g['A'])
PA0 = lam * (PA4 @ g['A'])

PA_2 = (PA0 + PA1 + PA2 + PA3 + PA4) / 5

lam = np.exp(-1j * 2 * 3*np.pi / 5)
PA1 = lam * (g['A'])
PA2 = lam * (PA1 @ g['A'])
PA3 = lam * (PA2 @ g['A'])
PA4 = lam * (PA3 @ g['A'])
PA0 = lam * (PA4 @ g['A'])

PA_3 = (PA0 + PA1 + PA2 + PA3 + PA4) / 5

lam = np.exp(-1j * 2 * 4*np.pi / 5)
PA1 = lam * (g['A'])
PA2 = lam * (PA1 @ g['A'])
PA3 = lam * (PA2 @ g['A'])
PA4 = lam * (PA3 @ g['A'])
PA0 = lam * (PA4 @ g['A'])

PA_4 = (PA0 + PA1 + PA2 + PA3 + PA4) / 5

# -- B projector
lam = np.exp(-1j * 2*np.pi / 4)
PB1 = lam * (g['B'])
PB2 = lam * (PB1 @ g['B'])
PB3 = lam * (PB2 @ g['B'])
PB0 = lam * (PB3 @ g['B'])

PB = (PB0 + PB1 + PB2 + PB3) / 3

# -- C projector
lam = np.exp(-1j * 2*np.pi / 2)
PC1 = lam * (g['AB'])
PC0 = lam * (PC1 @ g['AB'])

PC = (PC0 + PC1) / 2

HC = (id - 2*PC) + standard_coupling * Delta


n_x = 5

xs = np.linspace(0, 1, n_x)

path_AB = [[(1-x), x, 0] for x in xs]
path_BC = [[0, (1-x), x] for x in xs].pop(0)
path_CA = [[x, 0, (1-x)] for x in xs].pop(0)

path = path_AB + path_BC + path_CA

n_couplings = 100
n_energies = 200

couplings = np.linspace(0, 1, n_couplings)


op1 = (id - 2*PA_3)
op2 = Delta
# op2 = (id- 2*PA_2)
dosmat = np.zeros((n_couplings, n_energies))
for i in range(len(couplings)):
    H = couplings[i] * op1 + (1-couplings[i])*op2

    energy, dos = kpm.density_of_states(
        H, scale=5, n_moments=n_energies, n_energies=n_energies, n_random_states=10)

    dosmat[i] = dos

np.save("./energy.npy", energy)
np.save("./dos.npy", dosmat)


# dosmat = np.log10(1/dosmat+1)
plt.imshow(dosmat.transpose(), extent=(0, 1, np.amin(energy), np.amax(energy)), origin='lower', aspect='auto', cmap='Blues', interpolation='spline16',
           norm=colors.SymLogNorm(linthresh=0.1, linscale=0.5, vmin=0, vmax=np.amax(dosmat), base=10))
# plt.title("$H=(1-2P_A) + \epsilon \Delta$, $\mathrm{dim}~\mathcal{H} =$"+str(d))
# plt.title("$H=\epsilon (1-2P_A) + (1-\epsilon)\Delta$, $\mathrm{dim}~\mathcal{H} =$"+str(d))
# plt.colorbar(label=r"$\mathrm{dos}(E)$")
plt.colorbar()
ax = plt.gca()
# ax.set_xlabel(r'$\epsilon$')
# ax.set_ylabel(r'E')
ax.set_ylim((-4.2, 4.2))
plt.tight_layout()
plt.savefig('./dos.png', dpi=200)
# plt.show()


# PA = 1/5 * ( np.exp(-1j * np.pi/ 5) * g['A'].power(0) )
# PA = PA + PA.PA + PA.PA.PA + PA.PA.PA + PA.PA.PA.PA + PA.PA.PA.PA.PA

# energy, dos = kpm.density_of_states(Delta, scale=6.1, n_moments=100, n_energies=100, n_random_states=10)
# plt.plot(energy, dos)
# plt.show()
# plt.clf()
# spec = scipy.linalg.eigh(Delta.todense(), eigvals_only=True)

# plt.plot(spec)
# plt.show()
# H_fname = fname_prefix + ".hamiltonian"

# with open(H_fname, 'w') as H_file:
#     for i in range(d):
#         for j in H.rows[i]:
#             line = fortran_format(i) + " " + fortran_format(j) + " " + fortran_format(H[i,j]) + "\n"
#             H_file.write(line)
