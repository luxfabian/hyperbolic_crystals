"""
    ./py/tffg_group_extension.py

    Author: Fabian R. Lux
    Date:   9/11/2023
"""
import numpy as np

from scipy.sparse import diags, lil_matrix, kron

def get_extension(g=1, k=1):
    """
        Build and return the sigma matrices of the extension
    """

    s = np.exp(2*np.pi * 1j / k)

    # -- first sigma matrx_inv
    diagonal_0 = [ 1 for x in range(k) ]
    sigma_0 = lil_matrix(diags(diagonal_0, offsets=-1, shape=(k,k)))
    sigma_0[0,k-1] = 1

    # -- second sigma matrx_inv
    diagonal_1 = [ s**x for x in range(k) ]
    sigma_1 = lil_matrix(diags(diagonal_1, shape=(k,k)))

    # -- combine everything into a list
    sigma = np.array( [lil_matrix(np.eye(k)) for _ in range(4*g) ] )

    sigma[0] = sigma_0
    sigma[1] = sigma_1
    sigma[2*g] = sigma_0.transpose().conj()
    sigma[2*g+1] = sigma_1.transpose().conj()

    return sigma

def magnetic_representation(generators, k):
    """
        Represent the plain generators in the corresponding magnetic representation
    """

    # -- no. of generators
    n_g = len(generators)

    # -- genus
    g = n_g // 4

    sigma = get_extension(k)

    mag =  np.array( [ kron(sigma[j], generators[j]) for j in range(4*g) ] )

    return mag



if __name__ == '__main__':


    g = 2
    k = 5
    sigma = get_extension(g,k)

    # -- check commutation relations

    id = lil_matrix(np.eye(k))
    s = np.exp(2*np.pi * 1j / k) * id
    s_inv = np.exp(-2*np.pi * 1j / k) * id

    x = sigma[0]
    x_inv = sigma[2*g]
    y = sigma[0+1]
    y_inv = sigma[2*g+1]

    # s [x,y]
    relation = s.dot( x.dot(y.dot(x_inv.dot(y_inv)))) - id
    print("absmax(s [x,y] - 1) = ",np.amax( np.abs(np.array(relation)) ) )

    # [s,x]
    relation = s.dot(x.dot(s_inv.dot(x_inv))) - id
    print("absmax([s,x] - 1) = ",np.amax( np.abs(np.array(relation)) ) )

    # [s,y]
    relation = s.dot(y.dot(s_inv.dot(y_inv))) - id
    print("absmax([s,y] - 1) = ",np.amax( np.abs(np.array(relation)) ) )

    # s^k
    relation = np.linalg.matrix_power(s.toarray(),k) - id
    print("absmax(s^k - 1) = ",np.amax( np.abs(np.array(relation)) ) )