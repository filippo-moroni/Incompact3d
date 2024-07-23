#!---------------------------------------------------------!
#! With this script, we perform plotting of pre-multiplied !
#! 1D energy spectra of correlation functions.             !
#!---------------------------------------------------------!

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft

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

# Create folder to store plots
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
                                                                                                                                              
#!--------------------------------------------------------------------------------------!

#!--- Calculations ---!
          
# Inner quantities
sh_vel   = np.sqrt(nu * np.abs(mg_x[0]))  # shear velocity (based on streamwise mean gradient)  
delta_nu = nu / sh_vel                    # viscous length
y_plus   = y / delta_nu                   # y+

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
    
    # Friction Reynolds number
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

# Define wavenumber in spanwise direction (z)
kz = np.zeros(nz)
for i in range(len(kz)):
    kz[i] = (i+1)*(2*np.pi/Lz)

# FFT
Euuz = fft(Ruuz)


#!--- Plot section ---!

# Euuz
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = np.min(kz)
yliminf = np.max(kz)
xlimsup = np.min(Euuz)*1.2
ylimsup = np.max(Euuz)*1.2

# Euuz 
ax.scatter(kz, Euuz[c,:], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
# Axes labels
ax.set_xlabel(r'$k_z$',       fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$E_{uu}(z)$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Save and show the figure
save_and_show_plot('Euuz', snap_numb=snap_numb, add_string=add_string)









