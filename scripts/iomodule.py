"""
    ./scripts/iomodule.py

    Author: Fabian R. Lux
    Date:   2023-02-19

    Provides the routine get_access_point() which extracts information
    from the in and output files located at ${HYPERBOLIC_BUILD}/bin.

    This information can be used by other scripts to read the regular
    representation from file, that was previously generated in C++.
"""
import os
import re
import numpy as np

from scipy.sparse import lil_matrix

def get_access_point(silent=True):
    """
        Provides a dictionary access_point with access to the file
        locations of the regular representation and its dimensions.
        The keywords are:
            
            - 'p': p parameter of triangle group
            - 'q': q parameter of triangle group
            - 'N': finite approximation control parameter
            - 'd': dimension of the Hilbert space
            - 'fprefix': the prefix which should be used to ensure uniform file names
    """
    
    access_point = {}

    HYPERBOLIC_BUILD = os.getenv('HYPERBOLIC_BUILD')

    ###################################################################
    # Read group_specs.inp                                            #
    ###################################################################

    group_specs_fname = HYPERBOLIC_BUILD + "/bin/group_specs.inp"
    group_specs_file = open(group_specs_fname, 'r')
    group_specs_lines = group_specs_file.readlines()
    group_specs_lines = [ int(re.sub("[^0-9]", "", line)) for line in group_specs_lines]
    p,q,N = group_specs_lines

    access_point['p'] = p
    access_point['q'] = q
    access_point['N'] = N
    
    if not silent:
        print("Group specs: ", p, q, N) 

    group_specs_file.close()

    periodic_boundary = not (N<0)

    fname_prefix = HYPERBOLIC_BUILD + "/bin/{" + str(p) + "," + str(q) + "}_"
    if(periodic_boundary):
        fname_prefix += "modulo_"+str(N)
    else:
        fname_prefix += "open_"+str(-N) 

    access_point['fprefix'] = fname_prefix

    #######################################################################
    # Read spectral information                                           #
    #######################################################################

    spec_info_fname = fname_prefix + ".info"
    spec_info_file = open(spec_info_fname, 'r')
    spec_info_lines = spec_info_file.readlines()
    spec_info_lines = [ re.sub("[^0-9]", "", line) for line in spec_info_lines]
    d = int(spec_info_lines[0])

    access_point['d'] = d

    if not silent:
        print("Hilber space dimension: ", d) 
    
    spec_info_file.close()

    return access_point


    #######################################################################
    # Read regular represntation                                          #
    #######################################################################

#    generators = [ "A", "iA", "B", "iB" ]
#    H = lil_matrix((d,d))
#    for g in generators:
#        reg_fname = fname_prefix + "_" + g + ".reg" 
#        reg_data  = np.loadtxt(reg_fname, dtype=int, delimiter=" ")
#
#        for [i,j] in reg_data:
#            H[i,j] += 1.0
#
#
#    H_fname = fname_prefix + ".hamiltonian"
#
#    with open(H_fname, 'w') as H_file:
#      for i in range(d):
#        for j in H.rows[i]:
#          line = str(i) + " " + str(j) + " " + str(H[i,j]) + "\n"
#          H_file.write(line)
#
#    return access_point

if __name__=='__main__':

    access_point = get_access_point()

    print(access_point)














