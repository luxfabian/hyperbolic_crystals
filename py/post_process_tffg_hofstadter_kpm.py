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


from matplotlib.colors import SymLogNorm
from matplotlib.ticker import SymmetricalLogLocator

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


linthresh = 1.5*10**(0)  # The range within which the plot is linear
linscale = 1   # The factor by which data smaller than `linthresh` is scaled.

norm = SymLogNorm(linthresh=linthresh, linscale=linscale)
plt.imshow(hofstadter.transpose(), aspect='auto',norm=norm, extent=(0,1, Emin, Emax), interpolation='gaussian', origin = 'lower' ,resample=True,cmap='inferno')

ax = plt.gca()

# ax.set_yscale('symlog', linthresh=linthresh, linscale=linscale)

# The y-axis can sometimes have minor ticks in undesired places, this part ensures they're correctly placed
# ax.yaxis.set_minor_locator(SymmetricalLogLocator(linthresh=linthresh, base=10))


fs = 14
cbarfs = 12


ax.set_xlabel("$\phi$",fontsize=fs)
ax.set_ylabel("Energy",fontsize=fs)
ax.set_ylim( (-1,1))

cbar = plt.colorbar()
cbar.set_label(label=r'DOS (a.u.)', size=fs)


plt.title("$g=1$, $p^n=2^1$, $k=2003$, $n_{\mathrm{KPM}}=2048$, $n_{\mathrm{rand}}=10$",fontsize=11)
# plt.title("$\phi=1/{}$".format(k)+"")
plt.tight_layout()
plt.savefig("./out/hofstadter_kpm.png",dpi=1000)
plt.clf()