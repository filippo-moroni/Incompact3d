
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: This file stores Python subroutines for the pre-processing 
!              of Incompact3d simulations.
!              
!              Subroutines list:
!               - mem_and_cpuh.
!
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '../', '3-post_processing/python_common'))
sys.path.append(config_path)

#!-----------------------------------------------------------------------------!

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculation of:
!               - total number of points, n_tot;
!               - total number of snapshots, nsnap (one realization);
!               - total memory requirement [GB], in order to save velocity,
!                 pressure and 1 scalar field in double precision, accounting
!                 for all snapshots and all flow realizations.
!              Estimation of:
!               - total CPUh for the simulations (considering different
!                 flow realizations if present) assuming number of points 
!                 per CPU > 100'000.                  
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def mem_and_cpuh(nx,ny,nz,ifirst,ilast,itimescheme,ioutput,nrealiz):

    # Total number of points
    n_tot = nx*ny*nz

    # Calculate total number of snapshots (one realization)
    # Last '+1' is to count also the 1st time step (that is saved)
    nsnap = (ilast - ifirst + 1) // ioutput + 1

    # Total memory requirement [GB]: we are assuming to save velocity, pressure and 
    # one scalar field in double precision.
    mem_tot = nsnap * n_tot * 5 * 8.0 * (10**-9) * nrealiz
    mem_tot = round(mem_tot, 3)

    """
    Reference from ARIES runs of TTBL simulations:
    
    >>> Test 1:
     - ntot ~ 66 mln points;
     - CPUh ~ 13300;
     - points / nranks ~ 600'000;
     - ts_tot = 320'000;
     - results obtained with AB3 and implicit y-diffusion with CN.
    
    >>> Test 2:
     - ntot ~ 32 mln points;
     - CPUh ~ 1660;
     - points / nranks ~ 150'000;
     - ts_tot = 30'000;
     - results obtained with RK3 and implicit y-diffusion with CN.

    """
    
    # Adams-Bashforth 3
    if itimescheme == 3:
    
        n_tot_test  = 66 * (10**6)
        cpuh_test   = 13300
        ts_tot_test = 320000
    
    # Runge-Kutta 3
    elif itimescheme == 5:
    
        n_tot_test  = 32 * (10**6)
        cpuh_test   = 1660
        ts_tot_test = 30000
    
    
    # Ratio of CPUh of 1 flow realization by the product of number of points and number of total time-steps 
    performance_index = cpuh_test / (n_tot_test * ts_tot_test)  

    # Safety factor for estimation of total CPUh required
    # In ARIES, clock-frequency is 2.20 GHz.
    sf = 1.2  

    # Estimated CPUh (if we assume at least 100'000 points per core,
    # see 'Incompact3d user guide v1.1').
    cpuh = performance_index * n_tot * ilast * nrealiz * sf
    cpuh = round(cpuh, 1)
    
    return(n_tot, nsnap, mem_tot, cpuh)

#!-----------------------------------------------------------------------------!




