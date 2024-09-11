#!---------------------------------------------------------!
#! With this script, we perform plotting of pre-multiplied !
#! 1D energy spectra of correlation functions.             !
#!                                                         !
#! Inspired by 'spectra.py' by G. Boga.                    !
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
(itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, ioutput, iswitch_wo,  
 add_string, file1, filen, icrfile, nr, post_mean, post_vort, post_diss, post_corz, post_tke_eq
) = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!
    
#!--- Parameters ---!
uwall, nu, twd = set_flow_parameters(itype, re)
  
# Reading of grid points
y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)

# Read statistics data
(mean_u, mean_w, var_u, var_v, var_w, mean_uv, 
 vort_x, vort_y, vort_z, mg_tot, mg_x, mg_z,
 eps, Ruuz, Rvvz, Rwwz, Ruvz, Rssz,
 tke_conv, tke_turbt, tke_pstrain, tke_difft, tke_prod, tke_pseps,
 snap_numb) = read_data(itype, numscalar, post_mean, post_vort, post_diss, 
                        post_corz, post_tke_eq, ny, nz)
                                                                                                                                              
#!--------------------------------------------------------------------------------------!

#!--- Calculations ---!
          
# Inner quantities
sh_vel   = np.sqrt(nu * np.abs(mg_x[0]))  # shear velocity (based on streamwise mean gradient)  
delta_nu = nu / sh_vel                    # viscous length
y_plus   = y / delta_nu                   # y+
Ly_plus  = Ly / delta_nu                  # Ly+

delta_z  = Lz / nz

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
    re_tau = int(re_tau)
    
    # Print friction Reynolds number and boundary layer thickness
    print("Friction Reynolds number, re_tau = ", re_tau)
    print()
    print("Boundary layer thickness, delta_99 = ", bl_thick)
    print()
    print("Domain height in wall units, Ly+ = ", Ly_plus)
    print()

#!-------------------------------!

# Select between spectra and pre-multiplied spectra
i_premult = int(input("Would you like to plot pre-multiplied spectra? (0: no, 1: yes) "))
print()

#!--- y+ input for 1D spectra ---!

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

#!-------------------------------!

# Define angular wavenumber in spanwise direction (z)
# Since we have periodicity in z, the largest wavenumber is pi / delta_z 
# (smallest harmonic that can be represented (quarter of a sine, still a periodic function)
kz = np.zeros(nz)
for i in range(len(kz)):
    kz[i] = (i+1)*(2*np.pi/Lz)

# Define arrays 
Euuz = np.zeros((ny,nz))
Evvz = np.zeros((ny,nz))
Ewwz = np.zeros((ny,nz))
Euvz = np.zeros((ny,nz))

# Apply FFT (we need a periodic signal)
Euuz = np.real(fft(Ruuz))
Evvz = np.real(fft(Rvvz))
Ewwz = np.real(fft(Rwwz))
Euvz = np.real(fft(Ruvz))

# Half number of points in z to avoid periodicity effects
nzh = nz // 2

# Resize arrays and rescale in wall units 
kz     = kz[:nzh]*delta_nu
Euuz   = Euuz[:,:nzh] / sh_vel**2
Evvz   = Evvz[:,:nzh] / sh_vel**2
Ewwz   = Ewwz[:,:nzh] / sh_vel**2
Euvz   = Euvz[:,:nzh] / sh_vel**2

# Pre-multiply if asked
if i_premult == 1:
    Euuz = Euuz*kz
    Evvz = Evvz*kz
    Ewwz = Ewwz*kz
    Euvz = Euvz*kz

#!------------------------------------------------------------------------------------------------------------------------------!    
   
#!--- Plot section ---!

# List of variables and labels
variables = [('Euuz', 'E_{uu}^+(z)'), ('Evvz', 'E_{vv}^+(z)'), ('Ewwz', 'E_{ww}^+(z)'), ('Euvz', 'E_{uv}^+(z)')]

# Iterate over the variables to create plots
for var, ylabel in variables:

    # Subplot environment
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches, pp.yinches), linewidth=pp.tick_width, dpi=300)
    
    # Access the variable by name
    data = eval(var) 
    
    # Limits for axes
    xliminf = np.min(kz) * 0.8
    xlimsup = np.max(kz) * 1.2
    ylimsup = np.max(data[c, :]) * 1.2
    
    if i_premult == 1:
        yliminf = np.min(data) * 1.2
    else:
        yliminf = 0.000001

    # Plot
    ax.plot(kz, data[c, :], color='C0', linestyle='-', linewidth=pp.lw)

    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)
    
    # x-axis label
    ax.set_xlabel(r'$k_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    if i_premult == 1:
        
        # y-axis label
        ax.set_ylabel(r'$k_z^+{}$'.format(ylabel), fontsize=pp.fla, labelpad=pp.pad_axes_lab)
        
        # Save and show the figure
        save_and_show_plot(f'kz{var}', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, y_plus_in=y_plus_in)
    else:
    
        # y-axis label
        ax.set_ylabel(r'${}$'.format(ylabel), fontsize=pp.fla, labelpad=pp.pad_axes_lab)
        
        # Logarithmic y-axis (only in case of a standard spectrum Eij)
        ax.set_yscale('log')
        
        # Save and show the figure
        save_and_show_plot(var, snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, y_plus_in=y_plus_in)

#!------------------------------------------------------------------------------------------------------------------------------!







