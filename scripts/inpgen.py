"""
    ./scripts/inpgen.py
    
    Author: Fabian R. Lux
    Date:   2023-02-02
    
    Input generator which prepares the folder structure and basic input files
    for the execution of the C++ and Fortran programs.
"""

import os
import sys
import time

BUILD_DIRECTORY = os.environ['HYPERBOLIC_BUILD']
# OUT_DIRECTORY   = os.environ.get('HYPERBOLIC_OUT', BUILD_DIRECTORY+"\\out\\")

def mkdir(PATH):
    """
        Creates a directory at PATH if not already existent
    """
    if(os.path.isdir(PATH)):
        pass
    else:
        os.mkdir(PATH)

def set_unique_identifier():
    """
        Use the universal time, to create a unique identifier for the output folder.
    """

    # -- waiting for one second will ensure uniqueness of identifier
    time.sleep(1)

    OUT_DIR = os.getcwd() + 
    mkdir(OUT_DIR)

    if sys.platform in ['linux', 'macos']:

        HYPERBOLIC_DIR = OUT_DIR + '/' + str(round(time.time()))
        mkdir(HYPERBOLIC_DIR)

        # -- create .sh file to export HYPERBOLIC_DIR

        shellfile = HYPERBOLIC_DIR + '/hyperbolic_vars.sh'
        with open(shellfile,'w') as file:
            file.write("export HYPERBOLIC_DIR="+HYPERBOLIC_DIR)
            print(shellfile) 

    elif sys.platorm ['win32', 'cygwin']: 
        HYPERBOLIC_DIR = ".\\out\\" +str(round(time.time()))
        
    if(os.path.isdir(HYPERBOLIC_DIR)):
        #-- this should not happen
        pass
    else:
        os.mkdir(HYPERBOLIC_DIR)

    print("${HYPERBOLIC_DIR} =", HYPERBOLIC_DIR )
    

if __name__=="__main__":
    
    set_unique_identifier()
