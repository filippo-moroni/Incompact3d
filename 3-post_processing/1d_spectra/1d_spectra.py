#!---------------------------------------------------------!
#! With this script, we perform plotting of statistics     !
#! for TTBLs and channel flow simulations:                 !
#!                                                         !
#! - mean statistics (mean[u], var[u], etc.);              !
#! - mean total dissipation (to be done);                  !
#! - correlation coefficients for spanwise correlations.   !
#!                                                         !
#! Calculated and stored:                                  !
#!                                                         !
#! - non-dimensional grid spacings and domain dimensions.  !
#!---------------------------------------------------------!

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '..', 'python_common'))
sys.path.append(config_path)

# Import the plotting_params module
import plot_params as pp

# Import function to setting up, save and show plots 
from plot_subs import set_plot_settings, save_and_show_plot

# Import functions to read 'input.i3d', 'post.prm' files and statistics data 
from read_files import read_input_files, read_data

# Import function to setup flow parameters 
from set_flow_parameters import set_flow_parameters

#!--------------------------------------------------------------------------------------!

# Create folders to store later results (e.g. grid spacings and time scales files, plots)
os.makedirs('data_post', mode=0o777, exist_ok=True)
os.makedirs('plots',     mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
(itype, nx, ny, nz, istret, beta, 
 Lx, Ly, Lz, re, dt, numscalar, iswitch_wo, 
 add_string, file1, filen, icrfile, nr, 
 post_mean, post_vort, post_diss, post_corz, post_tke_eq) = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!
    
#!--- Parameters ---!
uwall, nu, twd = set_flow_parameters(itype, re)
  
# Reading of grid points
y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)

# Read statistics data
(mean_u, mean_w, var_u, var_v, var_w, mean_uv, 
 vort_x, vort_y, vort_z, mg_tot, mg_x, mg_z,
 eps, Ruuz, Rvvz, Rwwz, Ruvz, Rppz,
 tke_conv, tke_turbt, tke_pstrain, tke_difft, tke_prod, tke_pseps,
 snap_numb) = read_data(itype, numscalar, post_mean, post_vort, post_diss, 
                        post_corz, post_tke_eq, ny, nz)
                                                                                                                                   

# Valid only for TTBLs
if itype == 13:
    
    # Initialize the index
    j = 0
    
    # Calculate the index at which the BL thickness delta99 is and delta_99 itself
    while mean_u[j] > mean_u[0]*0.01: 
        
        # Boundary layer thickness delta_99
        bl_thick = y[j]
        
        # Increment the index
        j = j + 1
        
#!--------------------------------------------------------------------------------------!

#!--- Calculations ---!

# Mesh spacings
delta_x = Lx / nx
delta_z = Lz / nz
           
# Shear quantities
sh_vel   = np.sqrt(nu * np.abs(mg_x[0]))  # shear velocity (based on streamwise mean gradient)  
delta_nu = nu / sh_vel                    # viscous length
t_nu     = nu / (sh_vel ** 2)             # viscous time
  
# Rescaling variables through wall units
delta_x_plus = delta_x / delta_nu
delta_z_plus = delta_z / delta_nu
y_plus       = y       / delta_nu 

Lx_plus = Lx / delta_nu
Ly_plus = Ly / delta_nu 
Lz_plus = Lz / delta_nu

# Friction Reynolds number
if itype == 13:
    re_tau = sh_vel * bl_thick / nu
    
    # Print friction Reynolds number and boundary layer thickness
    print("Friction Reynolds number, re_tau = ", re_tau)
    print()
    print("Boundary layer thickness, delta_99 = ", bl_thick)
    print()
    print("Domain height in wall units, Ly+ = ", Ly_plus)
    print()


#!--- Calculations for correlations ---!

# Select the height at which correlations are plotted
y_plus_in = np.float64(input("Enter y+ value for correlations plotting: "))
print()

# Search for the index corresponding to the target y+ for correlations
c = 0 
for j in range(0, ny-1, 1):   
    if y_plus[j] < y_plus_in: c = c + 1

# Print the actual y+ value selected
print("Actual y+ value selected = ", y_plus[c])
print()

# Print the corresponding j-th index
print("Corresponding j-th index = ", c)
print()





