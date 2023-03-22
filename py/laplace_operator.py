"""
    ./py/laplace_operator.py

    Author: Fabian R. Lux
    Date:   2023-03-22

    Sets up the Laplace operator and saves it to file, such that it can be read in from Fortran.
"""

import iomodule
import numpy as np

from scipy.sparse import lil_matrix


def fortran_format(num):
    """
        Represent the number num (real or complex) as a Fortran readable string
        for E15.6 format.
    """

    fstring = ""   
    if isinstance(num, float):
        fstring = '{:15.6E}'.format(num)
    elif isinstance(num, complex):
        fstring = '{:15.6E}'.format(num.real) + " " + '{15.6E}'.format(num.imag)
    elif isinstance(num, int):
        fstring = '{:10}'.format(num)
    if "+" in fstring:
        fstring = fstring.replace("+", "") + " "
    return fstring.replace("E", "D")


def setup_laplacian():

    access_point = iomodule.get_access_point()



    d = access_point['d']
    fname_prefix = access_point['fprefix']

    generators = [ "A", "iA", "B", "iB" ]
    H = lil_matrix((d,d))
    for g in generators:
        reg_fname = fname_prefix + "_" + g + ".reg" 
        reg_data  = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

        for [i,j] in reg_data:
            H[i,j] += 0.25


    H_fname = fname_prefix + ".hamiltonian"

    with open(H_fname, 'w') as H_file:
        for i in range(d):
            for j in H.rows[i]:
                line = fortran_format(i) + " " + fortran_format(j) + " " + fortran_format(H[i,j]) + "\n"
                H_file.write(line)

if __name__=='__main__':
    setup_laplacian()