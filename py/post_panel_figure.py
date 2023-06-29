import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors

prl_figure_width =1*(3 + 3/8) #inch

aspect = 4
golden_ratio = 1.61803398875

# Some example data to display
x = np.linspace(-1, 1, 400)
y = np.sin(2*3.14*x ** 2)**2

fig, (ax1, ax2,ax3) = plt.subplots(1, 3,figsize=(aspect*prl_figure_width,prl_figure_width))
# fig.suptitle('Horizontally stacked subplots')

# ---------------------------------------------------------------------
#   IDS
# ---------------------------------------------------------------------

def integrated_density_of_states(spec, n_E):
    """
        Compute the integrated density of states from the spectral information
    """

    energy = np.linspace(-1,1,n_E)
    ids = np.zeros(n_E)

    for i in range(n_E):
        mu = energy[i]

        ids[i] = 0

        for E in spec:
            if E<mu:
                ids[i] += 1

    return energy, ids/len(spec)



files = ['../data/5_4_modulo_2.eig', '../data/5_4_modulo_3.eig', '../data/5_4_modulo_4.eig', '../data/5_4_modulo_5.eig', '../data/5_4_modulo_6.eig','../data/5_4_modulo_7.eig', '../data/5_4_modulo_8.eig']
mods =  [2,3,4,5,6,7,8]

sys_size = [160, 360,2560, 15000,57600,58800,81920]

files_subset = ['../data/5_4_modulo_2.eig','../data/5_4_modulo_4.eig','../data/5_4_modulo_8.eig']
mod_subset =  [2,4,8]
labels_subset = [r"$2,1,160=s,k, |G_k|$", r"$2,2,2560$" ,  r"$2,3,81920$" ]

labels = [r"$2,1,160=s,k, |G_k|$", r"$3,1,360$", r"$2,2,2560$" , r"$5,1,15000$" , r"$6,1,57600$" , r"$7,1,58800$", r"$2,3,81920$"  ]


ids_list = []



for i in range(len(files)):
    file = files[i]
    spec = np.load(file)
    energy, ids = integrated_density_of_states(spec,500)
    ids_list.append(ids)

    lw = 1
    if(i==len(files)-1):
        ax1.plot(energy, ids, label=labels[i], linewidth=lw,color='black')
        energy_last, ids_last = energy, ids
    else:
        ax1.plot(energy, ids, label=labels[i], linewidth=lw)

ax1.legend(fontsize=8)

ax1.text(-0.1, 1.15, "a)", transform=ax1.transAxes,
      fontsize=12, fontweight='bold', va='top', ha='right')

ax1.set_ylabel(r"Integrated density of states")
ax1.set_xlabel(r"Energy")
ax1.set_title(r"Periodic boundary condition")

def plot( pqN_tuple):

    # -- unpack tuple
    p, q, N = pqN_tuple

    # -- load spectrum
    spec = np.load("./out/{}_{}_open_{}.npy".format(p,q,N))

    energy, ids = integrated_density_of_states(spec,500)

    # -- plot
    ax2.plot(energy, ids  , '-',  label=r"{} sites".format(len(spec)) )


tuples  = [ (5,4,4),(5,4,6),(5,4,8),(5,4,10)]

for tuple in tuples:
    plot(tuple)

ax2.plot(energy_last,ids_last, '-',  label=r"$2,3,81920$",color='black') 

ax2.legend(fontsize=8,framealpha=1)

ax2.text(-0.1, 1.15, "b)", transform=ax2.transAxes,
      fontsize=12, fontweight='bold', va='top', ha='right')

ax2.set_ylabel(r"Integrated density of states")
ax2.set_xlabel(r"Energy")
ax2.set_title(r"Open boundary condition")
# ---------------------------------------------------------------------
#   fuchsian
# ---------------------------------------------------------------------

def plot( pn_tuple):

    # -- unpack tuple
    p, n = pn_tuple

    # -- load spectrum
    spec = np.load("./out/hyperbolic_chern_eig_{}_{}.npy".format(p,n))

    # -- plot
    ax3.plot(np.linspace(0,1,len(spec)),spec, 'o', markersize=1, label=r"${},{},{}$".format(p,n,int(len(spec)/4)) )

tuples  = [ (2,1), (2,2), (2,3), (2,4)]

for tuple in tuples:
    plot(tuple)

ax3.legend(fontsize=8,framealpha=1)


gap_1 = [-4.388824964976649,-2.345785512784003 ]
gap_2 = [-0.519871681082539,0.5264485364261997]
gap_3 = [2.352745485914674, 4.386258075803685]

id = np.ones(20)
# ca = ax3.gca()
# ax3.fill_between(np.linspace(0,1,20), gap_1[0]*id, gap_1[1]*id,color='lightgray')
ax3.fill_between(np.linspace(0,1,20), gap_2[0]*id, gap_2[1]*id,color='lightgray')
# ax3.fill_between(np.linspace(0,1,20), gap_3[0]*id, gap_3[1]*id,color='lightgray')

ax3.text(-0.1, 1.15, "c)", transform=ax3.transAxes,
      fontsize=12, fontweight='bold', va='top', ha='right')
ax3.set_xlabel(r"Relative index of eigenvalue")
ax3.set_ylabel(r"Energy")
ax3.set_title(r"Hyperbolic topological insulator candidate")

plt.tight_layout()
plt.savefig('./out/panel.png',dpi=300)
# plt.show()
plt.clf()