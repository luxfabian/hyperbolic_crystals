"""
    ./py/post_process_ABC_path.py

    Author: Fabian R. Lux
    Date:   2023-03-23

    Post process
"""

import numpy as np

fsize_major= 18
fsize_minor= 12

import matplotlib.pyplot as plt
import matplotlib.colors as colors

prl_figure_width =1*(3 + 3/8) #inch

aspect = 1.2
golden_ratio = 1.61803398875

energy = np.load("./out/ABCA_path_energy.npy")
dosmat = np.load("./out/ABCA_dos.npy")   

fig=plt.figure(figsize=(aspect*prl_figure_width,prl_figure_width))


im = plt.imshow(
        dosmat.transpose(), 
        extent=(0,1, np.amin(energy), np.amax(energy)), 
        origin='lower', aspect='auto',cmap='Blues', interpolation='spline16', 
        norm=colors.SymLogNorm(
            linthresh=0.1, linscale=0.6, vmin=0, vmax=np.amax(dosmat), base=10
            )
        ) 
    
cbar = plt.colorbar()

cbar.set_label(label='DOS', fontsize=fsize_minor)

cbar.ax.tick_params(labelsize=fsize_minor)

plt.title(r"Spectra of $\lambda_1 H_1 +  \lambda_2 H_2 + \lambda_3 H_3$",fontsize=fsize_major)
plt.grid(which='both')
ax=plt.gca()
ax.tick_params(labelsize=fsize_minor)
ax.set_ylabel("Spectrum",fontsize=fsize_major)
ax.set_ylim((-1.2,1.2))
#ax.set_xticks( [0, 1/6, 1/3, 1/3+1/6, 2/3, 2/3 + 1/6, 1] )
#ax.set_xticklabels( [r"$H_1$",r"$(H_1+H_2)/2$", r"$H_2$",r"$(H_2+H_3)/2$", r"$H_3$",r"$(H_3+H_1)/2$", r"$H_1$"])

ax.set_xticks( [0, 1/3,  2/3, 1] )

ax.set_xticks( [1/6,  1/3+1/6,  2/3 + 1/6] , minor=True)
ax.set_xticklabels( [r"$H_1$", r"$H_2$", r"$H_3$", r"$H_1$"],fontsize=16)
ax.set_xticklabels( [r"$(H_1+H_2)/2$", r"$(H_2+H_3)/2$",r"$(H_3+H_1)/2$"], minor=True,fontsize=8)
ax.set_xticklabels( ["","",""], minor=True,fontsize=8)



plt.tight_layout()
plt.savefig('./out/ABCA_path_dos.png',dpi=300)
# plt.show()
plt.clf()