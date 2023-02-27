"""
    ./py/hyperbolic_disk.py

    Author: Fabian R. Lux
    Date:   2023-02-27

    Basic routines which define how the triangle group acts on points in the hyperbolic plane.
"""
import numpy as np


def get_A(p, q):
    """
        Eq. (36) of PRB 105, 125118 (2022)
    """

    alpha = 2*np.pi / p

    A_mat = np.zeros((2, 2),dtype=complex)

    A_mat[0, 0] = np.exp(1j*alpha/2)
    A_mat[1, 1] = np.exp(-1j*alpha/2)

    return A_mat


def get_B(p, q):
    """
        Eq. (37) of PRB 105, 125118 (2022)
    """

    alpha = 2*np.pi / p
    beta = 2*np.pi / q
    r0 = get_r0(p, q)

    B_mat = np.zeros((2, 2),dtype=complex)

    B_mat[0, 0] = np.exp(1j*beta/2) - r0**2 * np.exp(-1j*beta/2)
    B_mat[1, 1] = np.exp(-1j*beta/2) - r0**2 * np.exp(+1j*beta/2)

    B_mat[0, 1] = r0 * (1-np.exp(1j*beta)) * np.exp(+1j*(alpha-beta)/2)
    B_mat[1, 0] = r0 * (1-np.exp(-1j*beta)) * np.exp(-1j*(alpha-beta)/2)

    B_mat = B_mat / (1-r0**2)

    return B_mat


def get_r0(p, q):
    """
        r0 as defined in Eq. (9) of PRB 105, 125118 (2022)
    """
    ca = np.cos(np.pi * (1/p + 1/q))
    cb = np.cos(np.pi * (1/p - 1/q))

    return np.sqrt(ca / cb)


def hyperbolic_distance(z1, z2):
    """
        Geodesic distance in the hyperbolic plane from z1 to z2
    """

    return np.arccosh(1 + 2 * abs(z1-z2)**2 / ((1 - abs(z1)**2) * (1 - abs(z2)**2)))


def moebius(sl2_mat, z):
    """
        The matrix sl2_mat acts from the left on the point z by a Moebius transformation.
        The image point is returned.
    """
    a = sl2_mat[0, 0]
    b = sl2_mat[0, 1]

    return (a*z + b) / (b.conj() * z + a.conj())

def parse_word(word, A, B, seed):
    """
        
    """

    z = seed
    for c in word[::-1]:
        if c=="A":
            z = moebius(A, z)
        elif c=="B":
            z = moebius(B, z)

    return z


if __name__=="__main__":

    p=5
    q=4

    A = get_A(p,q)
    B = get_B(p,q)

    seed = 0.4 + 0.1j

    print( parse_word("BBBB", A,B, seed) )