
"""
This file stores Python subroutines for the pre-processing 
of Incompact3d simulations. 
 
"""

#!-----------------------------------------------------------------------------!

# Libraries
import sys
import os
import numpy as np

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '../', '3-post_processing/python_common'))
sys.path.append(config_path)

# Import the plotting_params module
import plot_params as pp

#!-----------------------------------------------------------------------------!

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculation and estimation of:
!               - total number of points, n_tot;
!               - total number of snapshots, nsnap;
!               - total memory requirement [GB], in order to save velocity,
!                 pressure and 1 scalar field in double precision, accounting
!                 for all snapshots and all flow realizations;
!               - estimated total CPUh for the simulations
!                 (with points per CPU > 100'000).                  
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def mem_and_cpuh(nx,ny,nz,ifirst,ilast,ioutput,nrealiz):

    # Total number of points
    n_tot = nx*ny*nz

    # Calculate total number of snapshots
    # Last '+1' is to count also the 1st time step (that is saved)
    nsnap = (ilast - ifirst + 1) // ioutput + 1

    # Total memory requirement [GB]: we are assuming to save velocity, pressure and 
    # one scalar field in double precision.
    mem_tot = nsnap * n_tot * 5 * 8.0 * (10**-9) * nrealiz
    mem_tot = round(mem_tot, 3)

    """
    Reference from ARIES runs of TTBL simulations:
     - ntot ~ 66 mln points;
     - CPUh ~ 13300;
     - points / nranks ~ 600'000;
     - ts_tot = 320'000.
    """

    # Ratio of CPUh of 1 flow realization with the product of number of points and number of total time-steps 
    performance_index = 13000.0 / (66.0 * (10**6) * 320000)  

    # Safety factor for estimation of total CPUh required
    # In ARIES, clock-frequency is 2.20 GHz
    sf = 1.2  

    # Estimated CPUh (if we assume at least 100'000 points per core)
    cpuh = performance_index * n_tot * ilast * nrealiz * sf
    cpuh = round(cpuh, 1)
    
    return(n_tot, snap, mem_tot, cpuh)
    
 
     
