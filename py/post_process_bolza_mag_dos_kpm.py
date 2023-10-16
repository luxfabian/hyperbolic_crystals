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

def plot(param, label):

    N,n,k = param

    emesh = np.load("./out/bolza_"+str(N)+"_"+str(n)+"_"+str(k)+"_dos_kpm_emesh.npy")
    flux  = np.load("./out/bolza_"+str(N)+"_"+str(n)+"_"+str(k)+"_dos_kpm_flux.npy")
    dos = np.load("./out/bolza_"+str(N)+"_"+str(n)+"_"+str(k)+"_dos_kpm.npy")

    plt.plot(emesh,dos, label=label,linewidth=1)


# plot((4,3,10), label='$p^n=2^2, $\phi=1/1$')
# plot((3,1,3), label='$p^n=3^1$, $\phi=1/3$')
# plot((4,1,3), label='$p^n=2^2$')
# plot((3,1,3), label='$p^n=3^1$, $\phi=1/3$')
# plot((8,1,3), label='$p^n=2^3$')
# plot((4,1,1), label='$p^n=2^2$')
# plot((3,1,1), label='$p^n=3^1$,$\phi=0$')
# plot((8,1,1), label='$p^n=2^3$,$\phi=0$')
# plot((8,1,3), label='$p^n=2^3$,$\phi=1/3$')
plot((3,1,3), label='$p^n=2^3$,$\phi=1/3$')
plot((3,2,5), label='$p^n=2^3$,$\phi=2/5$')
# plot((8,1,1), label='$p^n=2^3$')
# plot((3,1,5), label='$p^n=3^1$, $\phi=1/3$')

plt.grid()

plt.title("$\phi=1/3$, $n_{\mathrm{KPM}}=512$, $n_{\mathrm{rand}}=10$")
# plt.legend(fontsize=12,framealpha=1,loc='lower right')

ax = plt.gca()

ax.legend(fontsize=12,framealpha=1)#,loc='upper left', bbox_to_anchor=(1, 1))

# ax.set_xlim((-1,-0.75))
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
plt.savefig("./out/bolza_dos_kpm.png",dpi=300)
plt.clf()