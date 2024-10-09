"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Small Python function to set up flow parameters and to read
!              y-coordinates of mesh elements: 'yp.dat'.
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
import numpy as np

def set_flow_parameters(re):

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

    # Kinematic viscosity (Reynolds is from the input file 'input.i3d')  
    nu = 1.0/re                       
        
    # Return to main program with parameters and y-coordinates
    return (uwall, nu, twd, y)
