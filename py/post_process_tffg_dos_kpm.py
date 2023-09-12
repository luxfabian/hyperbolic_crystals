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



def plot(param):

    N,k = param

    emesh = np.load("./out/tffg_"+str(N)+"_"+str(k)+"_dos_kpm_emesh.npy")
    flux  = np.load("./out/tffg_"+str(N)+"_"+str(k)+"_dos_kpm_flux.npy")
    dos = np.load("./out/tffg_"+str(N)+"_"+str(k)+"_dos_kpm.npy")

    plt.plot(-emesh,dos, label='$p^n={}$'.format(N), linewidth=0.5 + 0.1*N)


k = 5
plot((3,k))
plot((5,k))
plot((9,k))

plt.grid()

plt.title("$\phi=1/{}$".format(k)+", $n_{\mathrm{KPM}}=512$, $n_{\mathrm{rand}}=10$")
plt.legend(fontsize=12,framealpha=1)

ax = plt.gca()

# ax.set_xlim((-7,10))
# ax.set_ylim((0,0.01))

# plt.yscale("log")  

linthresh = 10**(-5)  # The range within which the plot is linear
linscale = 1   # The factor by which data smaller than `linthresh` is scaled.
ax.set_yscale('symlog', linthresh=linthresh, linscale=linscale)

# The y-axis can sometimes have minor ticks in undesired places, this part ensures they're correctly placed
ax.yaxis.set_minor_locator(SymmetricalLogLocator(linthresh=linthresh, base=10))


# The y-axis

ax.set_xlabel(r"Energy")
ax.set_ylabel(r"DOS [a.u.]")


plt.tight_layout()
plt.savefig("./out/tffg_dos_kpm.png",dpi=300)
plt.clf()