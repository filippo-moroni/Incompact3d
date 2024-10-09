
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: This file stores Python subroutines for the pre-processing 
!              of Incompact3d simulations.
!              
!              Subroutines list:
!               - mem_and_cpuh;
!               - plot_initial_vel_profile. 
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

# Import the plotting_params module
import plot_params as pp

# Import functions to setting up, save and show plots 
from plot_subs import set_plot_settings, save_and_show_plot

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

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Plot the initial velocity profile of the tanh initialization 
!              as Kozul et al. (2016). Calculate the shear layer initial 
!              thickness and the number of points contained in it.
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def plot_initial_vel_profile(ny,uwall,twd,theta_sl,yp,sh_vel_ic,delta_nu_ic):

    # Create folders to store later results
    os.makedirs('plots',                mode=0o777, exist_ok=True)
    os.makedirs('plots/pre_processing', mode=0o777, exist_ok=True)

    # Define the array
    Uo = np.zeros(ny)
        
    # Initial velocity profile (tanh) (Kozul et al. (2016))
    for j in range(0, ny):
        Uo[j] = uwall * (0.5 + 0.5 * (np.tanh((twd/2.0/theta_sl)*(1.0 - yp[j]/twd))))
    
    # Rescaling the initial velocity profile and the y-coordinates
    Uo = Uo / sh_vel_ic
    yp_ic = yp / delta_nu_ic

    # Calculate the thickness delta99^+ of the initial shear layer
    j = 0
    while j <= ny - 1 and Uo[j] > Uo[0]*0.01:
        sl_99_ic = yp_ic[j]
        j = j + 1
    	
    # Calculation of the number of mesh nodes in the initial shear layer
    npsl = 0      # number of points shear layer
    height = 0.0  # cumulative height in viscous unit (y+)
  
    for j in range(1, ny):
        if height + yp_ic[j] - yp_ic[j-1] <= sl_99_ic:
            npsl += 1 
        height += yp_ic[j] - yp_ic[j-1]

    """
    #!--- Plot section ---!

    # Plotting of the initial velocity profile in wall units, first 50 points
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
    
    # Description of .pdf file
    description = 'Initial streamwise velocity profile in wall units, as Kozul et al. (2016).'

    # Streamwise initial velocity profile
    ax.scatter(yp_ic, Uo, marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
    # Limits for axes
    xliminf = 0.0
    xlimsup = max(yp_ic)*1.2
    yliminf = min(Uo)*1.2
    ylimsup = max(Uo)*1.2
    
    # Axes labels
    ax.set_xlabel(r'$y^+$',   fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$U_o^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

    # Save and show the figure
    save_and_show_plot('init_umean', subfolder='pre_processing', description=description)
    
    """
    
    # Return calculated quantities
    return(sl_99_ic,npsl)



