"""
    ./py/tffg_.py

    Author: Fabian R. Lux
    Date:   9/11/2023
"""
import tffg_iomodule as iomodule
import numpy as np

import scipy

from tffg_group_extension import magnetic_representation

import matplotlib.pyplot as plt

n = 1

k = 1

access_point = iomodule.get_access_point()

g = access_point['g']
N = access_point['N']
d = access_point['d']

fname_prefix = access_point['fprefix']

# -- read generators

generators = []

for i in range(4*g):

    H = scipy.sparse.lil_matrix((d, d))

    reg_fname = fname_prefix + "_" + str(i) + ".reg"
    reg_data = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i, j] in reg_data:
        H[i, j] += 1

    generators.append(H)

generators = np.roll(generators, 2*g)
mag = magnetic_representation(generators, k,n)

# -- nearest neighbors
H_1_mag = []

# -- next-nearest neighbors
H_2_mag = []

for i in range(4*g): 
    H_1_mag.append(  (-mag[i])/(4*g)   )

    H_2_mag.append( (-mag[i].dot(mag[i]))/(4*g) )



dim = len(H_2_mag )

print("Dimension of parameter space: ", dim)

def get_gaps(spec):

    bandwidth = np.amax(spec) - np.amin(spec)

    gaps = ( (np.roll(spec,-1) - spec )[0:-1:1] ) / bandwidth

    return gaps

def cost_function(t):

    H_mag = scipy.sparse.csr_matrix(H_1_mag[0], dtype=complex)

    for i in range(len(H_1_mag)):

        H_mag += H_1_mag[i]

    for i in range(dim):

        H_mag += (t[2*i] + 1j*t[2*i+1])* H_2_mag[i]

    # -- hermitianize
    H_mag = (H_mag + H_mag.conj().transpose() ) / 2.0

    # -- diagonalize
    H_mag_dense = H_mag.todense()

    spec = np.sort(scipy.linalg.eigvalsh(H_mag_dense))#[1:-1:1]

    gaps = get_gaps(spec)

    cost = np.linalg.norm(gaps) 
    
    return cost, spec

 
def jacobian(t, h=1e-9):

    dim = len(H_2_mag)

    jac = np.zeros(2*dim,dtype=float)

    prev_cost, spec =  cost_function(t)

    for i in range(2*dim):
        
        dt = np.zeros(2*dim,dtype=float)
        dt[i] = h

        next_cost, _ = cost_function(t+dt)
        jac[i] = (next_cost - prev_cost ) / h

    return jac, spec

def step(t, eta=1):

    grad, spec = jacobian(t)

    err = np.linalg.norm(grad)

    return t + eta*grad, spec, err

if __name__=="__main__":

    # t = np.random.rand(2*dim)

    t = np.load("./out/tffg_gap_optimization_t.npy")

    index = 0
    eta = 2.5

    prev_err = 1e20

    prev_t = t
    while True:

        next_t, spec, next_err = step(prev_t, eta = eta)

        eta = 0.95*eta

        # if next_err > prev_err:
        #     # -- reject and adapt
        #     eta = 0.95*eta

        #     # print("Rejected")

        #     prev_t = next_t
        #     prev_err = next_err

        # else:
        prev_t = next_t
        prev_err = next_err
    
        plt.plot(spec,'o', markersize=1)

        # if index ==0:
        ax = plt.gca()
        ax.set_axisbelow(True)

        plt.grid()

        # ax.set_yscale("log")

        ax.set_xlabel("Eigenvalue")
        ax.set_ylabel("Energy")

        plt.tight_layout()

        print("Iteration no.", index)
        print("eta", eta)
        print("|| grad f || ", next_err)
        # print("t:" , t, "\n\n")

        index += 1
        plt.savefig("./out/tffg_gap_search_hist.png",dpi=300)
        plt.clf()

        np.save("./out/tffg_gap_optimization_t.npy", t)
        np.save("./out/tffg_gap_optimization_spec.npy", spec)

        if eta < 0.01:
            print("Max. iteration reached")
            break

  

        

    # res = scipy.optimize.minimize(cost_function,x0=t_guess)

    # jac = jacobian(t_guess)

    



    # H_1, H_2, H_1_mag, H_2_mag = operators()

    # spec, spec_mag = spectrum(H_1, H_2, H_1_mag, H_2_mag)

    # gaps = get_gaps(spec)
    # gaps_mag = get_gaps(spec_mag)

    # hist = np.histogram(gaps, bins=100, density=False)
    # hist_mag = np.histogram(gaps_mag, bins=100, density=False)

    # # plt.plot(np.sort(gaps), 'o')
   
    # plt.stairs(*hist, label="$\phi=0$")
    # plt.stairs(*hist_mag, label="$\phi=1/3$")

    # ax = plt.gca()
    # ax.set_axisbelow(True)

    # plt.legend()
    # plt.grid()

    # ax.set_yscale("log")
    # # ax.set_xlabel(r"$d(0, \sum_i  z_i  \| \psi_i \|^2 ) $")

    # ax.set_xlabel("$\Delta E / (E_{max} - E_{min})$")
    # ax.set_ylabel("Frequency")

    # plt.tight_layout()
    # plt.savefig("./out/tffg_gap_search_hist.png",dpi=300)
    # plt.clf()

    