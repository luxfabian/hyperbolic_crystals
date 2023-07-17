"""
    ./py/post_process_fuchsian_insulator.py

    Author: Fabian R. Lux
    Date:   6/27/2023

    We implement the model proposed in 

    https://doi.org/10.1038/s41467-023-36767-8

    for a hyperbolic Chern insulator with non-trivial second Chern number. We investigate,
    whether the band gaps of the system remain open in the thermodynamic limit.
"""

import numpy as np

fsize_major= 18
fsize_minor= 12

import matplotlib.pyplot as plt
import matplotlib.colors as colors

prl_figure_width =1*(3 + 3/8) #inch

aspect = 1.2
golden_ratio = 1.61803398875

# 5.134540750323415, -4.388824964976649
# 7.01940491591203, -2.345785512784003
# 8.661060802069857, -0.519871681082539
# 9.02587322121604, 0.5264485364261997
# 8.843467011642948, 2.352745485914674
# 9.086675291073739, 4.386258075803685

gap_1 = [-4.388824964976649,-2.345785512784003 ]
gap_2 = [-0.519871681082539,0.5264485364261997]
gap_3 = [2.352745485914674, 4.386258075803685]

# energy = np.load("./out/hyperbolic_chern_energy.npy")
# dos = np.load("./out/hyperbolic_chern_dos.npy")   

# fig=plt.figure(figsize=(aspect*prl_figure_width,prl_figure_width))


# plt.plot(energy,dos)

# plt.tight_layout()
# plt.savefig('./out/hyperbolic_chern_dos.png',dpi=300)
# # plt.show()
# plt.clf()

specs = {}
# specs[(2,1)] = np.sort( )
# spec_2_2 = np.sort( np.load("./out/hyperbolic_chern_eig_2_2.npy"))
# spec_3_1 = np.sort( np.load("./out/hyperbolic_chern_eig_3_1.npy"))
# spec_3_2 = np.sort( np.load("./out/hyperbolic_chern_eig_3_2.npy"))
# spec_2_3 = np.sort( np.load("./out/hyperbolic_chern_eig_2_3.npy"))
# spec_2_4 = np.sort( np.load("./out/hyperbolic_chern_eig_2_4.npy"))
# spec_5_1 = np.sort( np.load("./out/hyperbolic_chern_eig_5_1.npy"))
# spec_6_2 = np.sort( np.load("./out/hyperbolic_chern_eig_6_2.npy"))

# spec_test = np.sort( np.load("./out/hyperbolic_chern_eig_test.npy"))


fig=plt.figure(figsize=(aspect*prl_figure_width,prl_figure_width))


def plot( pn_tuple):

    # -- unpack tuple
    p, n = pn_tuple

    # -- load spectrum
    spec = np.load("./out/hyperbolic_chern_eig_{}_{}.npy".format(p,n))

    # -- plot
    plt.plot(np.linspace(0,1,len(spec)),spec, 'o', markersize=1.5, label=r"${},{},{}$".format(p,n,int(len(spec)/4)) )


# tuples  = [ (2,1), (3,1), (2,2), (6,1),  (5,1), (2,3), (7,1),(3,2), (11,1),(2,4), (13,1), (17,1)]

tuples  = [ (2,1), (3,1), (2,2), (6,1),  (5,1), (2,3),(3,2), (11,1),(2,4), (17,1)]

for tuple in tuples:
    plot(tuple)
# plt.plot(np.linspace(0,1,len(spec_test)), spec_test, 'o', markersize=1.5, label=r"test" )
# plt.plot(np.linspace(0,1,len(spec_2_2)), spec_2_2, 'o', markersize=0.5, label=r"$2,2,32$" )
# plt.plot(np.linspace(0,1,len(spec_2_2)), spec_2_2, 'o', markersize=0.5, label=r"$5,1,120$" )
# plt.plot(np.linspace(0,1,len(spec_2_3)),spec_2_3, 'o', markersize=0.5, label=r"$2,3,256$" )
# plt.plot(np.linspace(0,1,len(spec_3_2)),spec_3_2, 'o', markersize=0.5, label=r"$3,2,486$" )
# plt.plot(np.linspace(0,1,len(spec_2_4)),spec_2_4, 'o', markersize=0.5, label=r"$2,4,2048$" )
# plt.plot(np.linspace(0,1,len(spec_6_2)),spec_6_2, 'o', c='black', markersize=0.5, label=r"$6,2,7776$" )

plt.legend(fontsize=5,framealpha=1)
 
ax = plt.gca()

id = np.ones(20)

# ax.fill_between(np.linspace(0,1,20), gap_1[0]*id, gap_1[1]*id,color='lightgray')
ax.fill_between(np.linspace(0,1,20), gap_2[0]*id, gap_2[1]*id,color='lightgray')
# ax.fill_between(np.linspace(0,1,20), gap_3[0]*id, gap_3[1]*id,color='lightgray')
ax.set_ylabel(r"Energy")
ax.set_xlabel(r"Relative index of eigenvalue")

plt.tight_layout()
plt.savefig('./out/hyperbolic_spec.png',dpi=300)
# plt.show()
plt.clf()