"""
    ./py/fuchsian_iomodule.py

    Author: Fabian R. Lux
    Date:   2023-02-19

    Support for importing the Fuchsian group of genus two.

    Provides the routine get_access_point() which extracts information
    from the in and output files located at ${HYPERBOLIC_DIR}.

    This information can be used by other scripts to read the regular
    representation from file, that was previously generated in C++.
"""
import os
import re
import numpy as np


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

    HYPERBOLIC_DIR = os.getenv('HYPERBOLIC_DIR')

    ###################################################################
    # Read group_specs.inp                                            #
    ###################################################################

    group_specs_fname = HYPERBOLIC_DIR + "/fuchsian.inp"
    group_specs_file = open(group_specs_fname, 'r')
    group_specs_lines = group_specs_file.readlines()
    group_specs_lines = [int(re.sub('[^-\d]', "", line))
                         for line in group_specs_lines]
    p, n = group_specs_lines

    access_point['p'] = p
    access_point['n'] = n

    m = p**n

    if not silent:
        print("Group specs: ", p, n)

    group_specs_file.close()

   

    fname_prefix = HYPERBOLIC_DIR + "/fuchsian_" + str(m)
  
    access_point['fprefix'] = fname_prefix

    #######################################################################
    # Read spectral information                                           #
    #######################################################################

    spec_info_fname = fname_prefix + ".info"
    spec_info_file = open(spec_info_fname, 'r')
    spec_info_lines = spec_info_file.readlines()
    spec_info_lines = [re.sub("[^0-9]", "", line) for line in spec_info_lines]
    d = int(spec_info_lines[0])

    access_point['d'] = d

    if not silent:
        print("Hilber space dimension: ", d)

    spec_info_file.close()

    return access_point


if __name__ == '__main__':

    access_point = get_access_point()

    print(access_point)
