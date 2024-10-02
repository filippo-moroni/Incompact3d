"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Small Python function to set up flow parameters and to read
!              y-coordinates of mesh elements: 'yp.dat'.
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
import numpy as np

def set_flow_parameters(itype, re):

    # Try to read y-coordinates grid points from 'yp.dat'
    try:
    
        # Reading of y-coordinates grid points
        y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)
    
    # Handles the case where the file is not found (pre-processing)
    except OSError:  
        
        print("Warning: 'yp.dat' not found. Defaulting y to an empty array.")
        
        # Default to an empty array
        y = np.array([])  
        
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
    
    # Return to main program with parameters and y-coordinates
    return (uwall, nu, twd, y)
