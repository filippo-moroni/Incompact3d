#!-------------------------------------------!
#! Small Python function to set up           !
#! flow parameters                           !
#! (kinematic viscosity only at the moment). !
#!-------------------------------------------!

import numpy as np

def set_flow_parameters(itype, re):
 
    # Default parameters
    twd   = np.float64(1.0)                  # Trip wire diameter, D
    uwall = np.float64(1.0)                  # Wall velocity, Uwall 

    # TTBL
    if itype == 13:
    
        re    = np.float64(re)               # Trip Reynolds number
        nu    = 1.0/re                       # Kinematic viscosity

    # Channel (valid only for CFR at the moment)
    elif itype == 3:
        
        re_cent = np.float64(re)             # Centerline Reynolds number of the corresponding laminar Poiseuille flow
        re_tau  = 0.123*(re_cent**0.875)     # Estimated friction Reynolds number 
        nu      = 1.0/re_cent                # Kinematic viscosity
    
    # Return to main program with parameters
    return uwall, nu, twd
