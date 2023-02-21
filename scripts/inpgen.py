"""
    ./scripts/inpgen.py
    
    Author: Fabian R. Lux
    Date:   2023-02-02
    
    Input generator which prepares the folder structure and basic input files
    for the execution of the C++ and Fortran programs.
"""

import time
import os

BUILD_DIRECTORY =  os.environ['HYPERBOLIC_BUILD']
OUT_DIRECTORY = os.environ.get('HYPERBOLIC_OUT', BUILD_DIRECTORY+"\\out\\")

def set_unique_identifier():
    """
        Use the universal time, to create a unique identifier for the output folder.
    """
    time.sleep(1)
    
    HYPERBOLIC_DIR = OUT_DIRECTORY+str(round(time.time()))
    
    os.environ['HYPERBOLIC_DIR'] = HYPERBOLIC_DIR 
    
    if(os.path.isdir(HYPERBOLIC_DIR)):
        #-- this should not happen
        pass
    else:
        os.mkdir(HYPERBOLIC_DIR)
    
    print("${HYPERBOLIC_DIR} =", HYPERBOLIC_DIR )
    

if __name__=="__main__":
    
    set_unique_identifier()