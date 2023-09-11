"""
    ./py/post_process_tffg_mag.py

    Author: Fabian R. Lux
    Date:   9/11/2023
"""

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors

from dos import density_of_states

def process_spec(n,k):

    spec = np.load("./out/tffg_spec_"+str(n)+"_"+str(k)+".npy")


    for gamma in [0.2,0.15,0.1]:
        E_mesh, dos = density_of_states(spec, -7, 7, n_E=500, gamma=gamma)
        plt.plot(E_mesh, dos, '-', label='$\Gamma=${}'.format(gamma) )

    plt.grid()
    plt.title("$k={}$".format(k))
    plt.legend(fontsize=12,framealpha=1)
    
    ax = plt.gca()
    ax.set_xlabel(r"Energy")
    ax.set_ylabel(r"DOS [a.u.]")

    plt.tight_layout()
    plt.savefig("./out/tffg_spec_"+str(k)+".png",dpi=300)
    plt.clf()


# def compare_spec(n1,k1,n2,k2, gamma):

#     spec_1 = np.load("./out/tffg_spec_"+str(n1)+"_"+str(k1)+".npy")
#     E_mesh, dos_1 = density_of_states(spec_1, -7, 7, n_E=500, gamma=gamma)

#     spec_2 = np.load("./out/tffg_spec_"+str(n2)+"_"+str(k2)+".npy")
#     E_mesh, dos_2 = density_of_states(spec_2, -7, 7, n_E=500, gamma=gamma)

#     plt.plot(E_mesh, dos_1, '-', label='$k=${}'.format(k1) )
#     plt.plot(E_mesh, dos_2, '-', label='$k=${}'.format(k2) )

#     plt.grid()
#     plt.title("$\Gamma=${}".format(gamma))
#     plt.legend(fontsize=12,framealpha=1)
    
#     ax = plt.gca()
#     ax.set_xlabel(r"Energy")
#     ax.set_ylabel(r"DOS [a.u.]")

#     plt.tight_layout()
#     plt.savefig("./out/tffg_spec_"+str(k1)+"_"+str(k2)+".png",dpi=300)
#     plt.clf()

def compare_specs(nks, gamma,n_E=50):

    # plt.yscale("log")

    for nk in nks:
        N,n,k = nk

        

        spec = np.load("./out/tffg_spec_"+str(N)+"_"+str(n)+"_"+str(k)+".npy")

        print(n,k,len(spec))

        E_mesh, dos = density_of_states(spec, -7, 7, n_E=n_E, gamma=gamma)

        #plt.plot(E_mesh, dos, '-', linewidth=0.5 + 0.025*k, label='$\phi={}/{}$'.format(n,k) )
        plt.plot(E_mesh, dos, '-', linewidth=0.5 + 0.025*k, label='$p^n={}$'.format(N) )

    plt.grid()
    # plt.title("$g=2, p^n=2^3, \Gamma=${}".format(gamma))
    plt.title("$g=2, \phi=1/3, \Gamma=${}".format(gamma))
    plt.legend(fontsize=12,framealpha=1)
    
    ax = plt.gca()
    ax.set_xlabel(r"Energy")
    ax.set_ylabel(r"DOS [a.u.]")

    plt.tight_layout()
    plt.savefig("./out/tffg_compare_specs.png",dpi=300)
    plt.clf()
# process_spec(1)

# process_spec(3)

# nks = [(1,3),(2,6),(4,12),(8,24),(16,48)]

# tuples = [(2,1,3),(3,1,3),(4,1,3),(5,1,3),(6,1,3),(8,1,3)]

tuples = [(5,1,1),(5,1,3)]
compare_specs(tuples, 0.05,500)