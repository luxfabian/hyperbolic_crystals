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


count = 0

def plot(param):

    g,N,k = param

    emesh = np.load("./out/tffg_"+str(g)+"_"+str(N)+"_"+str(k)+"_dos_kpm_emesh.npy")
    flux  = np.load("./out/tffg_"+str(g)+"_"+str(N)+"_"+str(k)+"_dos_kpm_flux.npy")
    dos = np.load("./out/tffg_"+str(g)+"_"+str(N)+"_"+str(k)+"_dos_kpm.npy")

    plt.plot(-emesh,dos, label='$p^n={}$'.format(N), linewidth=1)



k = 3
# plot((1,2,k))
# plot((1,4,k))
# plot((1,8,k))
# plot((1,16,k))
# plot((1,32,k))
# plot((1,64,k))
plot((1,128,k))
plot((1,256,k))
plot((1,512,k))

plot((1,2,3001))

plt.grid()

plt.title("$g=1$, $\phi=1/{}$".format(k)+", $n_{\mathrm{KPM}}=2048$, $n_{\mathrm{rand}}=10$")
# plt.legend(fontsize=12,framealpha=1,loc='lower right')

ax = plt.gca()

ax.legend(fontsize=12,framealpha=1)#,loc='upper left', bbox_to_anchor=(1, 1))

ax.set_xlim((-1,1))
# ax.set_ylim((0,10))

# plt.yscale("log")  

# linthresh = 1.5*10**(-0)  # The range within which the plot is linear
# linscale = 1   # The factor by which data smaller than `linthresh` is scaled.
# ax.set_yscale('symlog', linthresh=linthresh, linscale=linscale)

# # The y-axis can sometimes have minor ticks in undesired places, this part ensures they're correctly placed
# ax.yaxis.set_minor_locator(SymmetricalLogLocator(linthresh=linthresh, base=10))


# The y-axis

fs = 14
ax.set_xlabel(r"Energy",fontsize=fs)
ax.set_ylabel(r"DOS (a.u.)",fontsize=fs)


plt.tight_layout()
plt.savefig("./out/tffg_dos_kpm.png",dpi=300)
plt.clf()