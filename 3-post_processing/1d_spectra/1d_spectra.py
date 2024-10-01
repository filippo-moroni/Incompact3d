
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: In this script we perform plotting of 1D energy spectra of 
!              correlation functions, with possibility to select between 
!              pre-multiplied or standard spectrum and to specify the 
!              y+ location.
!              Inspired by 'spectra.py' by G. Boga.
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# To be done: plotting of spectrum for scalar field correlations.

# Libraries
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '../../4-common', 'python_common'))
sys.path.append(config_path)

# Import the plotting_params module
import plot_params as pp

"""
Import functions to setting up, save and show plots and to obtain 
closest y+ location to what is chosen and its index
"""
from plot_subs import set_plot_settings, save_and_show_plot, y_plus_location

# Import functions to read 'input.i3d', 'post.prm' files and statistics data 
from read_files import read_input_files, read_data

# Import function to setup flow parameters 
from set_flow_parameters import set_flow_parameters

# Import function to calculate boundary layer thickness delta_99 for a TTBL
from ttbl_subs import calculate_ttbl_delta_99

#!--------------------------------------------------------------------------------------!

# Print to screen what the program does

print("!--- '1d_spectra.py' ---!")
print()
print(" Calculation and plotting of spectra of correlation coefficients")
print(" in spanwise direction at a certain y+ (selected from the user).")
print(" It is possible to plot pre-multiplied spectra if requested.")
print()

#!--------------------------------------------------------------------------------------!

# Create folder to store plots
os.makedirs('plots',         mode=0o777, exist_ok=True)
os.makedirs('plots/spectra', mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
(itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, ioutput, ioutput_cf, iswitch_wo,  
 add_string, file1, filen, icrfile, nr, post_mean, post_vort, post_diss, post_corz, post_tke_eq
) = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!
    
#!--- Parameters ---!
uwall, nu, twd = set_flow_parameters(itype, re)
  
# Reading of grid points
y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)

# Read statistics data
(mean_u, mean_w, var_u, var_v, var_w, mean_uv, 
 vort_x, vort_y, vort_z, mg_x, mg_z,
 eps, Ruuz, Rvvz, Rwwz, Ruvz, Rssz,
 tke_turbt, tke_presst, tke_difft, tke_prod, tke_pseps,
 snap_numb, i_switch_plot, ts, 
 sh_vel_x, sh_vel_tot) = read_data(itype, numscalar, post_mean, post_vort, post_diss, 
                                   post_corz, post_tke_eq, ny, nz)

#!--------------------------------------------------------------------------------------!

# Exit if the user tries to plot spectra using statistics from 'cf_monitoring'                                   
if i_switch_plot:
    print(">>> Spanwise correlation functions are present only in the dataset from 'post_incompact3d',")
    print(">>> we cannot calculate spectra and plot them.")
    print()    
    sys.exit(">>> Exiting from '1d_spectra.py'.")  
                                                                                                                                              
#!--------------------------------------------------------------------------------------!

#!--- Calculations ---!
          
# Inner quantities
delta_nu = nu / sh_vel_x    # viscous length
y_plus   = y  / delta_nu    # y+
Ly_plus  = Ly / delta_nu    # Ly+

# Spanwise mesh spacing
delta_z  = Lz / nz

# Valid only for TTBLs
if itype == 13:
    
    # Calculate BL thickness delta_99 for a TTBL and its related index
    (bl_thick, bl_thick_j) = calculate_ttbl_delta_99(mean_u, y)
    
    # Friction Reynolds number
    re_tau = sh_vel_x * bl_thick / nu
    re_tau = int(re_tau)

    # Print friction Reynolds number, boundary layer thickness and
    # domain height in viscous units
    print(">>> Friction Reynolds number, Re_tau = ", re_tau)
    print()
    print(">>> Boundary layer thickness, delta_99 = ", round(bl_thick,1))
    print()
    print(">>> Domain height in wall units, Ly+ = ", round(Ly_plus,1))
    print()

#!-------------------------------!

# Select between spectra and pre-multiplied spectra
i_premult = int(input(">>> Plot pre-multiplied spectra? (0: no, 1: yes) "))
print()

# Call external subroutine to determine the closest y+ location to what we want and its index 
(y_plus_index, y_plus_name) = y_plus_location(y_plus, ny)

#!-------------------------------!

# Define angular wavenumber in spanwise direction (z)
# Since we have periodicity in z, the largest wavenumber is pi / delta_z 
# (smallest harmonic that can be represented (quarter of a sine, still a periodic function).
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
Euuz   = Euuz[:,:nzh] / sh_vel_x**2
Evvz   = Evvz[:,:nzh] / sh_vel_x**2
Ewwz   = Ewwz[:,:nzh] / sh_vel_x**2
Euvz   = Euvz[:,:nzh] / sh_vel_x**2

# Pre-multiply if asked
if i_premult == 1:
    Euuz = Euuz*kz
    Evvz = Evvz*kz
    Ewwz = Ewwz*kz
    Euvz = Euvz*kz

#!------------------------------------------------------------------------------------------------------------------------------!    
   
#!--- Plot section ---!

# List of variables and labels
variables = [('Euuz', 'E_{uu}^+(z)', 'streamwise velocity auto-correlation in spanwise direction'       ), 
             ('Evvz', 'E_{vv}^+(z)', 'wall-normal velocity auto-correlation in spanwise direction'      ), 
             ('Ewwz', 'E_{ww}^+(z)', 'spanwise velocity auto-correlation in spanwise direction'         ), 
             ('Euvz', 'E_{uv}^+(z)', 'streamwise/wall-normal velocity correlation in spanwise direction')]

# Iterate over the variables to create plots
for var, ylabel, descr in variables:

    # Subplot environment
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches, pp.yinches), linewidth=pp.tick_width, dpi=300)
    
    # Access the variable by name
    data = eval(var) 
    
    # Limits for axes
    xliminf = np.min(kz) * 0.8
    xlimsup = np.max(kz) * 1.2
    ylimsup = np.max(data[y_plus_index, :]) * 1.2
    
    if i_premult == 1:
        yliminf = np.min(data) * 1.2
    else:
        yliminf = 0.000001

    # Plot
    ax.plot(kz, data[y_plus_index, :], color='C0', linestyle='-', linewidth=pp.lw)

    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)
    
    # x-axis label
    ax.set_xlabel(r'$k_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    if i_premult == 1:
        
        # y-axis label
        ax.set_ylabel(r'$k_z^+{}$'.format(ylabel), fontsize=pp.fla, labelpad=pp.pad_axes_lab)

        # Description of .pdf file
        description = 'Pre-multiplied spectrum of ' + descr + ' at y+ =' + f'{y_plus_name}' + '.'
        
        # Save and show the figure
        save_and_show_plot(f'kz{var}', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, y_plus_in=y_plus_name, subfolder='spectra', description=description)
    
    else:
    
        # y-axis label
        ax.set_ylabel(r'${}$'.format(ylabel), fontsize=pp.fla, labelpad=pp.pad_axes_lab)
        
        # Logarithmic y-axis (only in case of a standard spectrum Eij)
        ax.set_yscale('log')

        # Description of .pdf file
        description = 'Spectrum of ' + descr + ' at y+ =' + f'{y_plus_name}' + '.'
        
        # Save and show the figure
        save_and_show_plot(var, snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, y_plus_in=y_plus_name, subfolder='spectra', description=description)

#!------------------------------------------------------------------------------------------------------------------------------!







