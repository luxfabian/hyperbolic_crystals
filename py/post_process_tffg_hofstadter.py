"""
    ./py/post_process_tffg_mag.py

    Author: Fabian R. Lux
    Date:   9/11/2023
"""

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import RectBivariateSpline

from numba import njit

@njit
def fermi(E,mu,kbT):
    return 1.0 / ( np.exp( (E-mu) / kbT) + 1)

@njit
def dfermi(E,mu,kbT):
    return 1.0/(2*kbT + 2*kbT*np.cosh((E - mu)/kbT))
@njit
def dos(mu, kbT, spec):
       
    dens = 0.0 
    
    for level in spec:
        dens += dfermi(level,mu,kbT)
        
    return dens / ( len(spec) )

@njit
def ids(mu, spec):
       
    dens = 0.0 
    
    for level in spec:
        if level < mu:
            dens += 1
        
    return dens / ( len(spec) )




spec = np.load("./out/hofstadter.npy")
qs = np.load("./out/hofstadter_flux.npy")

# plt.imshow(hofstadter)

# plt.tight_layout()
# plt.savefig('./out/hofstadter.png',dpi=300)
# # plt.show()
# plt.clf()

nmu = 600
kbT = 0.02
n_q = len(qs)

Emins = np.zeros(n_q)
Emaxs  = np.zeros(n_q)
for i in range(n_q):
    Emins[i] = np.amin(spec[i])
    Emaxs[i] = np.amax(spec[i])
    
Emin = np.amin(Emins)
Emax = np.amax(Emaxs)

mus = np.linspace(Emin,Emax,nmu)


Q   = np.zeros(nmu*n_q)
DQ   = np.zeros(nmu*n_q)
DOS = np.zeros(nmu*n_q)
IDS = np.zeros(nmu*n_q)
MU  = np.zeros(nmu*n_q)

k = 0 
for i in range(nmu):
    for j in range(n_q):
        
        Q[k]   = qs[j]
        if j+1 < n_q:
            DQ[k]  = qs[j+1] - qs[j]
        IDS[k] = ids(mus[i], spec[j]) 
        DOS[k] = dos(mus[i], kbT, spec[j]) 
        MU[k]  = mus[i]
        
        k += 1

interpolation_order = 1
fids = RectBivariateSpline(mus, qs, np.reshape(IDS, (nmu,n_q)), kx=interpolation_order, ky=interpolation_order) 
fdos = RectBivariateSpline(mus, qs, np.reshape(DOS, (nmu,n_q)), kx=interpolation_order, ky=interpolation_order) 
fqdens = RectBivariateSpline(mus, qs, np.reshape(DQ, (nmu,n_q)), kx=interpolation_order, ky=interpolation_order) 


n_qi  = 500
n_mui = 500

qi = np.linspace(np.amin(qs),np.amax(qs),n_qi)

mui = np.linspace(Emin,Emax,n_mui)

Qi, MUi = np.meshgrid(qi,mui,indexing='ij')
DOSi = np.zeros( (n_qi, n_mui))
IDSi = np.zeros( (n_qi, n_mui))
QDENSi = np.zeros( (n_qi, n_mui))
for i in range(n_qi):
    for j in range(n_mui):
            
        DOSi[i,j] = np.abs(fdos( MUi[i,j], Qi[i,j]))
        IDSi[i,j] = np.abs(fids( MUi[i,j], Qi[i,j]))
        QDENSi[i,j] = np.abs(fqdens( MUi[i,j], Qi[i,j]))


fs = 14
cbarfs = 12

# DOSi2 = np.copy(DOSi)
# threshold = QDENSi > 0.25*maxdens
# DOSi2[threshold] = 0.0

# colmin = np.amin(DOSi2[DOSi2>0.001])
# colmax = np.amax(DOSi2[DOSi2>0.001])

plt.pcolormesh(Qi,MUi,DOSi, cmap='inferno')#,vmin=colmin, vmax=colmax);
ax = plt.gca()

ax.set_xlabel("$\phi=n/k$",fontsize=fs)
ax.set_ylabel("Energy",fontsize=fs)
ax.set_ylim( (Emin,Emax))

cbar = plt.colorbar()
cbar.set_label(label=r'DOS (a.u.)', size=fs)

# gaps = np.array( [ [0.2, 0],[0.8, 0.0], [0.3, 1.2], [0.7, 1.2], [0.3, -1.2], [0.7, -1.2] ] )

# gaps_class = [0,1,2,3,4,5]

# plt.scatter(gaps[:,0],gaps[:,1],c=gaps_class, cmap='tab20c',s=10)



plt.title("$g=2$, $p^n=2^1$, $k=200$",fontsize=fs)

plt.tight_layout()
plt.savefig('./out/hofstadter_dos.png',dpi=300)
# plt.show()
plt.clf()

# -- ids

plt.pcolormesh(Qi, IDSi, MUi, cmap='RdBu')#,vmin=-0.5,vmax=7)
ax = plt.gca()
ax.set_xlabel("$\phi=n/k$",fontsize=fs)
ax.set_ylabel("IDS",fontsize=fs)
# ax.set_ylim(1,1.8)
cbar = plt.colorbar()
cbar.set_label(label="Energy", size=fs)


# # tups = [  [0,1],[0,2],[0,3],[0,4],[0,5],[1,-1],[1,-2],[1,-3],[1,-4],[1,-5],[1,1],[1,2],[1,3],[1,4],[1,5] ]
# tups = []

# for n in np.arange(-5,6,1):
#     for k in np.arange(-5,6,1):
#         tups.append([n,k])

# for tup in tups:
#     n , k = tup
#     spec2 = ( n/4 + k*qs /4  ) 
#     plt.plot(qs,spec2, color='black', linewidth=.5, linestyle='-')

# ax.set_ylim((0,1))
# plt.scatter(gaps_q, gaps_ids,c=gaps_class, cmap='Wistia', s=20, zorder=10)

# ax.annotate('1,0,1', xy=(0.35, 1.08), rotation=0, fontsize=10)
# ax.annotate('1,0,2', xy=(0.35, 1.20), rotation=0, fontsize=10)
# ax.annotate('1,0,3', xy=(0.35, 1.33), rotation=0, fontsize=10)
# ax.annotate('2,-2,2', xy=(0.32, 1.48), rotation=0, fontsize=10)
plt.title("$g=2$, $p^n=2^1$, $k=200$",fontsize=fs)
plt.tight_layout()
plt.savefig('./out/hofstadter_ids.png',dpi=300)
plt.clf() 