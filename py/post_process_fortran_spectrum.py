"""
    ./py/fortran_spec.py

    Author: Fabian R. Lux
    Date:   2023-03-22

    Load the spectrum which was computed in Fortran and post-process the result. 
    This will generate Fig. 2 of the manuscript
"""

import iomodule
import numpy as np
import matplotlib.pyplot as plt

prl_figure_width = 1.5*(3 + 3/8) #inch
golden_ratio = 1.61803398875

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

def process_spectrum():

    access_point = iomodule.get_access_point()

    # eig_fname = access_point['fprefix'] + '.eig'

    eig_fname = '../data/5_4_modulo_5.eig'

    spec = np.load(eig_fname)

    print(np.amin(spec))

    energy, ids = integrated_density_of_states(spec,200)

    plt.plot(energy, ids)
    plt.show()

def multi_ids_plot():

    files = ['../data/5_4_modulo_2.eig', '../data/5_4_modulo_3.eig', '../data/5_4_modulo_4.eig', '../data/5_4_modulo_5.eig', '../data/5_4_modulo_6.eig','../data/5_4_modulo_7.eig', '../data/5_4_modulo_8.eig']
    mods =  [2,3,4,5,6,7,8]

    sys_size = [160, 360,2560, 15000,57600,58800,81920]

    files_subset = ['../data/5_4_modulo_2.eig','../data/5_4_modulo_4.eig','../data/5_4_modulo_8.eig']
    mod_subset =  [2,4,8]
    labels_subset = [r"$2,1,160=s,k, |G_k|$", r"$2,2,2560$" ,  r"$2,3,81920$" ]

    labels = [r"$2,1,160=s,k, |G_k|$", r"$3,1,360$", r"$2,2,2560$" , r"$5,1,15000$" , r"$6,1,57600$" , r"$7,1,58800$", r"$2,3,81920$"  ]
    ids_list = []

    fig = plt.figure(figsize=(prl_figure_width, prl_figure_width/golden_ratio),dpi=300,layout='constrained')

    for i in range(len(files)):
        file = files[i]
        spec = np.load(file)
        energy, ids = integrated_density_of_states(spec,500)
        ids_list.append(ids)

        lw = 1
        if(i==len(files)-1):
            plt.plot(energy, ids, label=labels[i], linewidth=lw,color='black')
        else:
            plt.plot(energy, ids, label=labels[i], linewidth=lw)

        

    plt.legend(fontsize=8)

    ax = plt.gca()
    ax.set_ylabel(r"Integrated density of states")
    ax.set_xlabel(r"Energy")

    # center zoom
    # ax.set_xlim((-0.1,0.1))
    # ax.set_ylim((0.4,0.6))

    plt.savefig("./out/ids_convergence.png",dpi=300)
    #plt.show()
    plt.close()

    ref = ids_list[-1]

    Err = []

    dE = energy[1]-energy[0]

    norm = np.trapz(ref, dx=dE) 

    for i in range(len(files)-1):

        err= np.trapz( np.abs( ids_list[i] - ref )**2 , dx=dE) / 2

        Err.append(err)

    s = 2.6
    fig = plt.figure(figsize=(1.1*prl_figure_width/s, prl_figure_width/s), dpi=300,layout='constrained')

    plt.plot(sys_size[:-1],Err,'-')
    plt.plot(sys_size[:-1],Err,'o',color='black')
    

    ax = plt.gca()

    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top') 
    
    #ax.set_xticks((2,3,4,5,6,7))
    #ax.set_xticklabels(( "(2,1)" ,"(3,1)","(2,2)","(5,1)","(6,1)","(7,1)"),fontsize=7)
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.set_ylim( (0.35*10**-5, 10**-3))
    ax.set_ylabel(r"MSE")
    ax.set_xlabel(r"$|G_k|$")
    #ax.set_aspect(1.6)
    plt.grid()
    #plt.tight_layout()
    plt.savefig("./out/ids_convergence_mse.png",dpi=300, transparent=True)
    #plt.show()

def ground_state_convergence():

    files = ['../data/5_4_modulo_2.eig','../data/5_4_modulo_3.eig', '../data/5_4_modulo_4.eig', '../data/5_4_modulo_5.eig', '../data/5_4_modulo_6.eig', '../data/5_4_modulo_7.eig', '../data/5_4_modulo_8.eig']

    E = []
    for file in files:
        spec = np.load(file)

        E.append(np.amin(spec))


    print(E)

    
    mods = [2,3,4,5,6,7,8]
    plt.plot(mods,E,'o')
    plt.title(r"Ground state energy of $\Delta ~\mathrm{mod}~ m$")
    plt.grid()
    ax  = plt.gca()

    ax.set_xlabel(r"$m$")
    ax.set_ylabel(r"$E$")

    ax.set_xticks( (1,2,3,4,5,6,7,8))
    plt.tight_layout()
    plt.show()

        

if __name__=='__main__':

    # process_spectrum()

    # ground_state_convergence()

    multi_ids_plot()