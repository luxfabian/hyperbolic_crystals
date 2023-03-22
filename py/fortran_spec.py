"""
    ./py/fortran_spec.py

    Author: Fabian R. Lux
    Date:   2023-03-22

    Load the spectrum which was computed in Fortran and post-process the result
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

    files = ['../data/5_4_modulo_2.eig', '../data/5_4_modulo_3.eig', '../data/5_4_modulo_4.eig', '../data/5_4_modulo_5.eig', '../data/5_4_modulo_6.eig']
    mods =  [2,3,4,5,6]

    labels = [r"$2,1,160=s,k, |G_k|$", r"$3,1,360$", r"$2,2,2560$" , r"$5,1,15000$" , r"$6,1,57600$"  ]
    ids_list = []

    fig = plt.figure(figsize=(prl_figure_width, prl_figure_width/golden_ratio),dpi=300,layout='constrained')

    for i in range(len(files)):
        file = files[i]
        spec = np.load(file)
        energy, ids = integrated_density_of_states(spec,200)
        ids_list.append(ids)

        lw = 1.5
        if(i==len(files)-1):
            lw = 1.5*lw

        plt.plot(energy, ids, label=labels[i], linewidth=lw)

    plt.legend(fontsize=8)

    ax = plt.gca()
    ax.set_ylabel(r"Integrated density of states")
    ax.set_xlabel(r"Energy")

    plt.savefig("ids_convergence.png",dpi=300)
    #plt.show()
    plt.close()

    ref = ids_list[-1]

    Err = []

    dE = energy[1]-energy[0]

    norm = np.trapz(ref, dx=dE)

    for i in range(len(files)-1):

        err= np.trapz( np.abs( ids_list[i] - ref )**2 , dx=dE) / 2

        Err.append(err)

    fig = plt.figure(figsize=(1*prl_figure_width/3, prl_figure_width/3), dpi=300,layout='constrained')

    plt.plot(mods[:-1],Err,'-')
    plt.plot(mods[:-1],Err,'o',color='black')
    

    ax = plt.gca()

    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top') 
    
    ax.set_xticks((2,3,4,5))
    ax.set_xticklabels(( "(2,1)" ,"(3,1)","(2,2)","(5,1)"),fontsize=7)
    ax.set_yscale('log')
    ax.set_ylim( (0.35*10**-5, 10**-3))
    ax.set_ylabel(r"MSE")
    ax.set_xlabel(r"$(s,k)$")
    ax.set_aspect(1.35)
    plt.grid()
    #plt.tight_layout()
    plt.savefig("ids_convergence_mse.png",dpi=300, transparent=True)
    #plt.show()

def ground_state_convergence():

    files = ['../data/5_4_modulo_2.eig','../data/5_4_modulo_3.eig', '../data/5_4_modulo_4.eig', '../data/5_4_modulo_5.eig', '../data/5_4_modulo_6.eig']

    E = []
    for file in files:
        spec = np.load(file)

        E.append(np.amin(spec))


    print(E)

    
    mods = [2,3,4,5,6]
    plt.plot(mods,E,'o')
    plt.title(r"Ground state energy of $\Delta ~\mathrm{mod}~ m$")
    plt.grid()
    ax  = plt.gca()

    ax.set_xlabel(r"$m$")
    ax.set_ylabel(r"$E$")

    ax.set_xticks( (1,2,3,4,5,6))
    plt.tight_layout()
    plt.show()

        

if __name__=='__main__':

    # process_spectrum()

    # ground_state_convergence()

    multi_ids_plot()