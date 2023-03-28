"""
    ./fuchsian_group/kpm.py

    Author: Fabian R. Lux
    Date:   28/12/2022

    Implementation of the kernel polynomial method. 
    Code is adapted and boiled down from:

        https://github.com/joselado/kpmpy

    The implementation here relies on vectorization and just-in-time
    compilation instead of Fortran.

    See DOI: 10.1103/RevModPhys.78.275 for a detailed description of 
    the method itself.
"""

import numpy as np
import scipy

from numba import njit


@njit
def jackson_kernel(mus):
    """
        Jackson damping of high-frequency oscillations
    """

    n = len(mus)
    pn = np.pi/(n+1)

    tan = np.tan(pn)

    i = np.arange(n)
    cos = np.cos(pn*i)
    sin = np.sin(pn*i)

    fac = ((n-i+1)*cos+sin/tan)/(n+1)

    return mus*fac


@njit
def lorentz_kernel(mus):
    """
        Lorentz damping of high-frequency oscillations
    """

    n = len(mus)

    i = 3*(1-np.arange(n)/n)
    fac = np.sinh(i) / np.sinh(3)

    return mus*fac


@njit
def chebyshev_reconstruction(mus, xs, kernel="jackson"):
    """
        Reconstructs a function at the abcissas stored in `xs`, using 
        its Chebyshev moments stored in `mus`. A damping kernel is
        applied to the moments before the reconstruction. The Chebyshev
        polynomials are defined iterativeley.

        mus:    Chebyshev moments
        xs:     Abcissas where function is to be reconstructed
        kernel: Damping kernel
    """

    # -- storage of Chebyshev polynomials (T0(x)=1)
    tm = np.zeros(xs.shape) + 1.
    t = xs.copy()

    # -- Kernel selection
    if kernel == "jackson":
        mus = jackson_kernel(mus)
    elif kernel == "lorentz":
        mus = lorentz_kernel(mus)
    else:
        raise Exception

    # -- Chebyshev recursion
    ys = np.zeros(xs.shape) + mus[0]
    for i in range(1, len(mus)):
        mu = mus[i]
        ys += 2.*mu*t  # add contribution
        tp = 2.*xs*t - tm  # chebychev recursion relation
        tm = t + 0.
        t = 0. + tp  # next iteration

    ys = ys/np.sqrt(1.-xs*xs)/np.pi
    return ys


def chebyshev_moments(v, m, n_moments=100):
    """
        n-th order Chebychev expansion of <v|m|v>.

        m is rescaled such that Spec(m) C [-1,1]
    """
    v = v.todense()
    v = np.array(v).reshape(v.shape[0])

    # -- moments
    mus = np.zeros(2*n_moments, dtype=complex)

    am = v.copy()  # zero vector
    a = m*v  # vector number 1
    bk = (np.conjugate(v).T)@v  # scalar product
    bk1 = (np.conjugate(a).T)@v  # scalar product

    mus[0] = bk.copy()  # mu0
    mus[1] = bk1.copy()  # mu1
    for i in range(1, n_moments):
        ap = 2*m@a - am  # recursion relation
        bk = (np.conjugate(a).T)@a  # scalar product
        bk1 = (np.conjugate(ap).T)@a  # scalar product
        mus[2*i] = 2.*bk
        mus[2*i+1] = 2.*bk1
        am = a.copy()  # new variables
        a = ap.copy()  # new variables
    mu0 = mus[0]  # first
    mu1 = mus[1]  # second
    for i in range(1, n_moments):
        mus[2*i] += -mu0
        mus[2*i+1] += -mu1
    return mus


@njit
def random_state(d):
    """
        Create a normalized random state in the Hilbert space of
        dimension d.

        d: Hilbert space dimension
    """

    # -- entries are uniformly distributed in [-0.5,0.5)
    state = np.random.rand(d) - 0.5

    # -- normalize
    state = state/np.linalg.norm(state)

    return state


def random_trace(m, n_random_states=20, n_moments=200):
    """ 
        Approximate the trace of m by sampling over random states.

        m:                  csc_matrix
        n_random_states:    no. of random states
        n:                  no. of Chebyshev moments
    """

    # -- Hilbert space dimension
    d = m.shape[0]  # length of the matrix

    # -- Chebyshev momentschebyshev_moments
    mus = np.zeros(2*n_moments, dtype=complex)

    for _ in range(n_random_states):

        # -- create random state
        v = scipy.sparse.csc_matrix(random_state(d)).transpose()

        # -- evaluate moments
        mus += chebyshev_moments(v, m, n_moments=n_moments)

    return mus/n_random_states


def density_of_states(H, scale=1., n_moments=10, n_energies=500,
                      n_random_states=20, epsilon=0.01, kernel="jackson"):
    """
        KPM calculation of the density of states

        H:                      Hamiltonian (csc_matrix)
        scale:                  rescaling such that Spec(H/scale) C [-1,1]
        n_moments:              no. of moments
        n_energies:             no. of energy values to sample the spectrum
        kernel:                 damping kernel 
        n_random_states:        no. of random states for trace evaluation
        epsilon:                spectral regularization parameter
    """

    # -- rescaled Hamiltonian such that Spec(m) C [-1,1]
    m = H/scale

    # -- energy values for DOS reconstruction
    E = np.linspace(-1, 1, n_energies) * (1-epsilon)

    # -- obtaine moments
    moments = random_trace(
        m, n_random_states=n_random_states, n_moments=n_moments)

    # -- reconstruct dos
    dos = chebyshev_reconstruction(moments, E, kernel=kernel).real

    return (scale*E, dos/scale)


if __name__ == '__main__':
    pass
