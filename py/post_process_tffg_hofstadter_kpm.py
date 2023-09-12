"""
    ./py/post_process_tffg_mag.py

    Author: Fabian R. Lux
    Date:   9/11/2023
"""

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from scipy.interpolate import RectBivariateSpline

from numba import njit


emesh = np.load("./out/hofstadter_kpm_emesh.npy")
flux = np.load("./out/hofstadter_kpm_flux.npy")
hofstadter = np.load("./out/hofstadter_kpm.npy")

Emin = np.amin(emesh)
Emax = np.amax(emesh)

phimin = np.amin(flux)
phimax = np.amax(flux)

# plt.imshow(hofstadter.transpose(), aspect='auto', extent=(phimin,phimax, Emin, Emax), interpolation='bilinear')

# plt.imshow(hofstadter.transpose(), aspect='auto',norm=LogNorm(), extent=(phimin,phimax, Emin, Emax), interpolation='bilinear')

plt.imshow(hofstadter.transpose(), aspect='auto',norm=LogNorm(), extent=(phimin,phimax, Emin, Emax), interpolation='gaussian', origin = 'upper' )

fs = 14
cbarfs = 12
ax = plt.gca()

ax.set_xlabel("$\phi=n/k$",fontsize=fs)
ax.set_ylabel("Energy",fontsize=fs)
ax.set_ylim( (Emin,Emax))

cbar = plt.colorbar()
cbar.set_label(label=r'DOS (a.u.)', size=fs)


plt.title("$g=2$, $p^n=3^1$, $k=101$",fontsize=fs)

plt.tight_layout()
plt.savefig("./out/hofstadter_kpm.png",dpi=300)
plt.clf()