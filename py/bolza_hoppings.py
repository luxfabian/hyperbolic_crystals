"""
    ./py/bolza_hoppings.py

    Author: Fabian R. Lux
    Date:   10/14/2023

    Implements the hopping matrices of the proper {8,3}-tesselation 
    of the Bolza surface.
                    
    These matrices can be used to decorate the generators of the Fuchsian
    translation group \Gamma_2
"""
import numpy as np

from scipy.sparse import lil_matrix


def bolza_hoppings():
    """
        Construct the hopping matrices in sparse format.

        Returns three objects TA0, TB0, TA1:

        TA0: intra-unit cell A-hoppings
        TB0: intra-unit cell B-hoppings
        TA1: inter-unit cell A-hoppings
    """

    # -- auxiliary basis lookup

    d = 48
    labels = np.zeros((d, 3), dtype=int)
    states = np.zeros((8, 2, 3), dtype=int)
    i = 0
    for k in range(8):
        for l in range(2):
            for m in range(3):
                states[k, l, m] = i
                labels[i] = [k, l, m]
                i += 1

    # -- intracell B-cycles

    TB0 = lil_matrix((d, d))

    for j in range(d):

        k, l, m = labels[j]

        i = states[k, l, (m - 1) % 3]

        TB0[i, j] += 1


    # -- intracell A-cycles

    TA0 = lil_matrix((d, d))

    for mu in range(8):

        i = states[(mu - 1) % 8, 0, 0]
        j = states[mu, 0, 0]

        TA0[i, j] += 1

        i = states[(mu + 1) % 8, 1, 0]
        j = states[(mu + 1) % 8, 0, 1]

        TA0[i, j] += 1

        i = states[(mu + 1) % 8, 0, 1]
        j = states[mu, 0, 2]

        TA0[i, j] += 1

        i = states[mu, 0, 2]
        j = states[mu, 1, 2]

        TA0[i, j] += 1

    # -- intercell A-cycles

    TA1 = [lil_matrix((d, d)) for _ in range(8)]
    for mu in range(8):

        i = states[ (mu + 4) % 8, 1, 2]
        j = states[ (mu+1) % 8, 1, 0]

        TA1[mu][i, j] += 1

        i = states[(mu + 5) % 8, 1, 1]
        j = states[mu, 1, 1]

        TA1[mu][i, j] += 1

    return TA0, TB0, TA1


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    TA0, TB0, TA1 = bolza_hoppings()

    plt.imshow(TA0.toarray(),cmap="Blues")
    plt.colorbar()
    plt.title(r"$T_{A,0}$")
    plt.savefig("TA0.png",dpi=300)
    plt.clf()

    plt.imshow(TB0.toarray(),cmap="Blues")
    plt.colorbar()
    plt.title(r"$T_{B,0}$")
    plt.savefig("TB0.png",dpi=300)
    plt.clf()

    for mu in range(8):
        plt.imshow(TA1[mu].toarray(),cmap="Blues")
        plt.colorbar()
        plt.title("$T_{A,"+str(mu+1)+"}$")
        plt.savefig("TA{}.png".format(mu+1),dpi=300)
        plt.clf()

    pass
