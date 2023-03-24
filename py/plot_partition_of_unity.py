import numpy as np
from hyperbolic_disk import *
import matplotlib.pyplot as plt

n_phi = 301
n_r = 100

phis = np.linspace(0, 2*np.pi, n_phi)
rs = np.linspace(0,0.99999,n_r)

chis = np.zeros((3,n_r, n_phi))

phi = np.pi/5
for i in range(n_r):
    for j in range(n_phi):
        z = rs[i] * np.exp(1j*phis[j])
        for region in range(3):
            chis[region,i,j] = smooth_phase(z,1.0,phi=phi,region=region)

fig, axs = plt.subplots(1,3, subplot_kw=dict(projection="polar"),figsize=(2,0.8),dpi=300,layout='constrained')
deg = 180 / np.pi

p1,p2,p3 = junction_angles(phi)

titles = [r"$\chi_1(z)$", r"$\chi_2(z)$" , r"$\chi_3(z)$"]
for region in range(3):
    ax = axs[region]
   # ax.set_xticklabels([r"$\Upsilon_1^\infty$", r"$\Upsilon_2^\infty$", r"$\Upsilon_3^\infty$"],fontsize=5)
    ax.set_xticklabels(["","",""],fontsize=5)
    ax.set_yticklabels([])
    ax.set_thetagrids([p1*deg, p2*deg, p3*deg])
    ax.set_rticks([])
    ax.set_title(titles[region],fontsize=10)
    ctf = ax.contourf(phis,rs,chis[region], levels=400, cmap='viridis')

    
plt.tight_layout(w_pad=1)
sk = 1
fig.colorbar(ctf, ax=axs.ravel().tolist(), orientation='vertical', ticks=[0,0.5,1], shrink=sk, aspect=5*sk, pad=0.1)
plt.savefig("Y_junction_1.0.png", dpi=300) 
plt.show()