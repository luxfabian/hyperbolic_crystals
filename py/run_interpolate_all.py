"""
    ./py/run_interpolate_all.py

    Author: Fabian R. Lux
    Date:   2023-02-28

    Interpolates between various topological Hamilton operators and 
    calculates the spectrum along the path using the KPM method
"""

# -- options -----------------------------

import matplotlib.colors as colors
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
import numpy as np
import scipy
import kpm
import iomodule
n_couplings = 200
n_energies = 200
n_moments = n_energies
n_random_states = 10

# ----------------------------------------


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
Delta = (g['A'] + g['iA'] + g['B'] + g['iB'])/4

id = scipy.sparse.identity(d)

# -- A projectors
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

PB_1 = (PB0 + PB1 + PB2 + PB3) / 4

lam = np.exp(-1j * 2 * 2*np.pi / 4)
PB1 = lam * (g['B'])
PB2 = lam * (PB1 @ g['B'])
PB3 = lam * (PB2 @ g['B'])
PB0 = lam * (PB3 @ g['B'])

PB_2 = (PB0 + PB1 + PB2 + PB3) / 4

lam = np.exp(-1j * 3 * 2*np.pi / 4)
PB1 = lam * (g['B'])
PB2 = lam * (PB1 @ g['B'])
PB3 = lam * (PB2 @ g['B'])
PB0 = lam * (PB3 @ g['B'])

PB_3 = (PB0 + PB1 + PB2 + PB3) / 4

# -- C projector
lam = np.exp(-1j * 2*np.pi / 2)
PC1 = lam * (g['AB'])
PC0 = lam * (PC1 @ g['AB'])

PC_1 = (PC0 + PC1) / 2


projectors = [PA_1, PA_2, PA_3, PA_4, PB_1, PB_2, PB_3, PC_1]
labels = ["PA_1", "PA_2", "PA_3", "PA_4", "PB_1", "PB_2", "PB_3", "PC_1"]


couplings = np.linspace(0, 1, n_couplings)

for i in range(len(labels)):

    # path = "/mnt/c/Users/flux/Dropbox/research/projects/2022_03_nyc/Emil+Fabian/RegularTessellations/Simulations/spectral_interpolations"
    path = "."
    fname = path+"/out/"+labels[i] + "_" + str(d) + "_"

    op1 = id-2*projectors[i]
    op2 = Delta

    dosmat = np.zeros((n_couplings, n_energies))

    for i in range(len(couplings)):
        H = couplings[i] * op1 + (1-couplings[i])*op2

        energy, dos = kpm.density_of_states(
            H, scale=1.2, n_moments=n_moments, n_energies=n_energies, n_random_states=n_random_states)

        dosmat[i] = dos

    np.save(fname+"energy.npy", energy)
    np.save(fname+"dos.npy", dosmat)

    fig = plt.figure()
    ax = plt.gca()
    ax.set_ylim((-1.2, 1.2))
    ax.set_xticks([0, 0.25, 0.50, 0.75, 1.00])

    plt.imshow(
        dosmat.transpose(),
        extent=(0, 1, -1.2, 1.2),
        origin='lower', aspect='auto', cmap='Blues', interpolation='spline16',
        norm=colors.SymLogNorm(
            linthresh=0.1, linscale=0.6, vmin=0, vmax=np.amax(dosmat), base=10
        )
    )

    plt.colorbar()

    plt.tight_layout()
    plt.savefig(fname+'dos.png', dpi=200)
    plt.clf()
