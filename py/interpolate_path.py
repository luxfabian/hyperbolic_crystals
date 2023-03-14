"""
    ./scripts/projector_interpolation.py

    Author: Fabian R. Lux
    Date:   2023-02-28

    Interpolates between to projection operators and calculates the spectrum along the path
"""

# -- options -----------------------------

n_x = 100
n_energies = 200
n_moments = 80
n_random_states = 10

# ----------------------------------------

import iomodule
import kpm
import scipy
import numpy as np
from scipy.sparse import lil_matrix

import matplotlib.pyplot as plt
import matplotlib.colors as colors

plt.rcParams.update({'font.size': 14})

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
Delta = (g['A'] + g['iA'] + g['B'] + g['iB'])/4

id = scipy.sparse.identity(d)

# -- A projectors
lam = np.exp(-1j * 2*np.pi/ 5)
PA1 = lam *  ( g['A'] )
PA2 = lam *  ( PA1 @ g['A'] )
PA3 = lam *  ( PA2 @ g['A'] )
PA4 = lam *  ( PA3 @ g['A'] )
PA0 = lam *  ( PA4 @ g['A'] )

PA_1 = ( PA0 + PA1 + PA2 + PA3 + PA4) / 5

lam = np.exp(-1j *2* 2*np.pi/ 5)
PA1 = lam *  ( g['A'] )
PA2 = lam *  ( PA1 @ g['A'] )
PA3 = lam *  ( PA2 @ g['A'] )
PA4 = lam *  ( PA3 @ g['A'] )
PA0 = lam *  ( PA4 @ g['A'] )

PA_2 = ( PA0 + PA1 + PA2 + PA3 + PA4) / 5

lam = np.exp(-1j *2* 3*np.pi/ 5)
PA1 = lam *  ( g['A'] )
PA2 = lam *  ( PA1 @ g['A'] )
PA3 = lam *  ( PA2 @ g['A'] )
PA4 = lam *  ( PA3 @ g['A'] )
PA0 = lam *  ( PA4 @ g['A'] )

PA_3 = ( PA0 + PA1 + PA2 + PA3 + PA4) / 5

lam = np.exp(-1j *2* 4*np.pi/ 5)
PA1 = lam *  ( g['A'] )
PA2 = lam *  ( PA1 @ g['A'] )
PA3 = lam *  ( PA2 @ g['A'] )
PA4 = lam *  ( PA3 @ g['A'] )
PA0 = lam *  ( PA4 @ g['A'] )

PA_4 = ( PA0 + PA1 + PA2 + PA3 + PA4) / 5

# -- B projector
lam = np.exp(-1j * 2*np.pi/ 4)
PB1 = lam *  ( g['B'] )
PB2 = lam *  ( PB1 @ g['B'] )
PB3 = lam *  ( PB2 @ g['B'] )
PB0 = lam *  ( PB3 @ g['B'] )

PB_1 = ( PB0 + PB1 + PB2 + PB3 ) / 4

lam = np.exp(-1j * 2* 2*np.pi/ 4)
PB1 = lam *  ( g['B'] )
PB2 = lam *  ( PB1 @ g['B'] )
PB3 = lam *  ( PB2 @ g['B'] )
PB0 = lam *  ( PB3 @ g['B'] )

PB_2 = ( PB0 + PB1 + PB2 + PB3 ) / 4

lam = np.exp(-1j * 3* 2*np.pi/ 4)
PB1 = lam *  ( g['B'] )
PB2 = lam *  ( PB1 @ g['B'] )
PB3 = lam *  ( PB2 @ g['B'] )
PB0 = lam *  ( PB3 @ g['B'] )

PB_3 = ( PB0 + PB1 + PB2 + PB3 ) / 4

# -- C projector
lam = np.exp(-1j * 2*np.pi/ 2)
PC1 = lam *  ( g['AB'] )
PC0 = lam *  ( PC1 @ g['AB'] )

PC_1 = ( PC0 + PC1 ) / 2

epsilon = 0.8
HA = epsilon*(id - 2 * PA_1) + (1-epsilon) * Delta
HB = epsilon*(id - 2 * PB_1) + (1-epsilon) * Delta
HC = epsilon*(id - 2 * PC_1) + (1-epsilon) * Delta

xs = np.linspace(0,1,n_x)

path_AB = [ [(1-x),x,0] for x in xs ]
path_BC = [ [0, (1-x),x] for x in xs ]
path_CA = [ [x, 0, (1-x)] for x in xs ]

path_BC.pop(0)
path_CA.pop(0)

path = path_AB + path_BC + path_CA

n_pts = len(path)

dosmat = np.zeros( (n_pts, n_energies) )
for i in range(n_pts):
    H = path[i][0] * HA + path[i][1] * HB + path[i][2] * HC
    energy, dos = kpm.density_of_states(H, scale=1.2, n_moments=n_moments, n_energies=n_energies, n_random_states=n_random_states)

    dosmat[i] = dos

np.save("./out/ABCA_path_energy.npy", energy)
np.save("./out/ABCA_dos.npy", dosmat)       

plt.imshow(
        dosmat.transpose(), 
        extent=(0,1, np.amin(energy), np.amax(energy)), 
        origin='lower', aspect='auto',cmap='Blues', interpolation='spline16', 
        norm=colors.SymLogNorm(
            linthresh=0.1, linscale=0.6, vmin=0, vmax=np.amax(dosmat), base=10
            )
        ) 
    
plt.colorbar(label='DOS')
ax=plt.gca()
ax.set_ylabel("Energy")
ax.set_ylim((-1.2,1.2))
ax.set_xticks( [0, 1/3, 2/3, 1] )
ax.set_xticklabels( [r"$H_1$", r"$H_2$", r"$H_3$", r"$H_1$"])
plt.tight_layout()
plt.savefig('./out/ABCA_path_dos.png',dpi=200)
plt.show()
plt.clf()
