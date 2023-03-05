"""
    ./py/hyperbolic_disk.py

    Author: Fabian R. Lux
    Date:   2023-02-27

    Basic routines which define how the triangle group acts on points in the hyperbolic plane.
"""
import numpy as np

def ahlfors(z1,z2):
    """
        Defines the Ahlfors bracket
    """
    return (1-abs(z1)**2) * (1- abs(z2)**2)  + abs(z1-z2)**2

def midpoint(z1,z2):
    """
        Calculate the geodesic midpoint between z1 and z2 
        in the hyperbolic disk.
    """
    
    nom  = z1 * (1-abs(z2)**2) + z2 * (1-abs(z1)**2)
    denom = 1 - abs(z1)**2 * abs(z2)**2 + ahlfors(z1,z2)* np.sqrt( (1-abs(z1)**2)*(1-abs(z2)**2))

    return nom/denom

def phase(z, phi):
    """
	Determines the phase of z
    """

    arg = np.arctan2(z.imag, z.real) % (2*np.pi)

    p1 = phi
    p2 = phi + 2*np.pi/3
    p3 = phi + 4*np.pi/3

    phase = 2
    if p1 <= arg and arg < p2:
        phase = 0
    if p2 <= arg and arg < p3:
        phase = 1

    return phase 

def sigmoid(x):
    return 1/(1+np.exp(x))

def smooth_phase(z,l,phi,region):
    chiz = chi(z,l,phi)
  
    return chiz[region] 

def get_d(z,alpha):
    dist = hyperbolic_distance(0,z)
    if abs(alpha) >= np.pi/2:
        return(dist)
    else:
        return np.arcsinh( abs(np.sin(alpha)) * np.sinh(dist))

def chi(z,l,phi):
    
    # -- determine in which phase I am in
    this_phase = phase(z, phi)
    
    # -- initialize phase boundaries
    dphi = 2*np.pi/3
    Y = np.zeros(3, dtype=complex)
    Y[0] = np.exp(1j*phi)
    Y[1] = np.exp(1j*dphi) * Y[0]
    Y[2] = np.exp(1j*dphi) * Y[1]

    # -- calculate angles
    def arg(z):
        return np.arctan2(z.imag,z.real)

    ds = np.zeros(3)
    if this_phase == 0:
        ds[1] = get_d(z, arg(Y[1]*z.conj()))
        ds[2] = get_d(z, arg(Y[0]*z.conj()))
        ds[0] = - min(ds[1],ds[2]) 
    if this_phase == 1:
        ds[0] = get_d(z, arg(Y[1]*z.conj()))
        ds[2] = get_d(z, arg(Y[2]*z.conj()))
        ds[1] = - min(ds[2],ds[0]) 
    if this_phase == 2:
        ds[0] = get_d(z, arg(Y[0]*z.conj()))
        ds[1] = get_d(z, arg(Y[2]*z.conj()))
        ds[2] = - min(ds[0],ds[1]) 
        
    sigma  = [ sigmoid(d/l) for d in ds]
  
    norm = sigma[0] + sigma[1] + sigma[2]

    return sigma/norm

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
    import matplotlib.pyplot as plt

    n_phi = 301
    n_r = 100

    phis = np.linspace(0, 2*np.pi, n_phi)
    rs = np.linspace(0,0.99999,n_r)

    chis = np.zeros((3,n_r, n_phi))
 
    phi = np.pi/5
    for i in range(n_r):
        for j in range(n_phi):
            z = rs[i] * np.exp(1j*phis[j])
            for region in range(3):
                chis[region,i,j] = smooth_phase(z,0.5,phi=phi,region=region)

    fig, axs = plt.subplots(1,3, subplot_kw=dict(projection="polar"))
    deg = 180 / np.pi
  
    titles = [r"$\chi_1(z)$", r"$\chi_2(z)$" , r"$\chi_3(z)$"]
    for region in range(3):
        ax = axs[region]
       	ax.set_xticklabels([])
       	ax.set_yticklabels([])
        ax.set_thetagrids([phi*deg, phi*deg+120, phi*deg+240])
        ax.set_rticks([])
        ax.set_title(titles[region])
        ctf = ax.contourf(phis,rs,chis[region], levels=400, cmap='viridis')
    plt.tight_layout()
    
    sk = 0.4
    fig.colorbar(ctf, ax=axs.ravel().tolist(), orientation='vertical', ticks=[0,0.5,1], shrink=sk, aspect=20*sk)
    plt.savefig("Y_junction.png", dpi=300) 
    plt.show()
