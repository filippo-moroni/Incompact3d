#!--------------------------------------------------------!
#! With this script, we perform calculations and plotting !
#! 1D pre-multiplied spectra.of statistics                !
#!                                                        !
#! Inspired by 'spectra.py' by G. Boga                    !
#!--------------------------------------------------------!

# to be completed

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

# Import function to setup flow parameters (kinematic viscosity only at the moment)
from set_flow_parameters import set_flow_parameters

#!--------------------------------------------------------------------------------------!

# Create folders to store later results (e.g. grid spacings and time scales files, plots)
os.makedirs('data_post', mode=0o777, exist_ok=True)
os.makedirs('plots',     mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
itype, nx, ny, nz, Lx, Ly, Lz, re, numscalar, iswitch_wo, file1, filen, icrfile, nr, add_string = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!

# Settings for contourf and colormap

cmap_name     = "Greys"
field_label   = r'$k_z^+ E_{uu}^+$'
xlabel        = r'$k_z$'
pad_cbar_lab  = -8
size_cbar     = '10%'
    
#!--- Parameters ---!
uwall, nu = set_flow_parameters(itype, re)
 
# Reading of grid points
y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)
y = y[:ny]

# Read statistics data
(mean_u, mean_w, var_u, var_v, mean_uv, 
 vort_x, vort_y, vort_z, mg_tot, mg_x, mg_z,
 eps, Ruuz, Rvvz, Rwwz, Ruvz, Rppz,
 snap_numb) = read_data(itype, numscalar)
 
# Calculate y+
sh_vel   = np.sqrt(nu * np.abs(mg_x[0]))  # shear velocity (based on streamwise mean gradient)  
delta_nu = nu / sh_vel                    # viscous length 
y_plus   = y / delta_nu                                                                                   

# Define variables
kz     = np.zeros(nz)
kzEuuz = np.zeros((ny,nz))

# Create wavenumber array in z-dir.
for k in range(len(kz)):
    kz[k] = (k+1)*(2.0*np.pi/Lz)
    kz[k] = kz[k] / delta_nu

# Multiply correlation functions by the wavenumber and perform FFT 
for j in range(0, ny, 1):
    for k in range(0, nz-1, 1):
        Ruuz[j,k] = Ruuz[j,k] / sh_vel
        Ruuz[j,k] = Ruuz[j,k] * kz[k]
    kzEuuz[j,:] = np.float64(fft(Ruuz[j,:]))

#!--- Plot 1D section ---!

# Limits for axes (used in 'set_plot_settings')
xliminf = 0.1
xlimsup = 1000.0
yliminf = 0.0
ylimsup = y_plus[-1]

# Extent of the image (dimensions of the domain)
extent = [xliminf, xlimsup, yliminf, ylimsup]

#!--------------------------------------------------------------------------------------!

#!--- Mesh section and iso-levels ---!

# Values of iso-levels        
lvls = np.linspace(np.min(kzEuuz), np.max(kzEuuz), pp.nlvl)

field_ticks = [np.min(kzEuuz),np.max(kzEuuz)]

X, Y = np.meshgrid(kz, y_plus)

#!--------------------------------------------------------------------------------------!

# 1D pre-multiplied spectrum of longitudinal velocity correlations in spanwise direction

# Subplots environment
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
    
# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

# Functions to locate the colorbar
divider = make_axes_locatable(ax)
    
cax = divider.append_axes('right', size=size_cbar, pad=0.05)
        
# Imshow function (unexpectedly it adjusts well the aspect ratio of the plotted image with contourf)
#im = ax.imshow(kzEuuz, cmap=cmap_name, extent=extent, origin='upper')
                
# Plotting with filled contours    
C = ax.contourf(X, Y, kzEuuz, lvls, cmap=cmap_name, extend='neither')
    
# Colorbar
cbar = fig.colorbar(C, cax=cax, orientation='vertical', ticks=field_ticks)
       
# Colorbar ticks 
cbar.ax.tick_params(axis='y', labelsize=pp.fla2, length=pp.lmajt, width=pp.tick_width) 
     
# Colorbar label (use pp.pad_cbar_lab to use the default value for padding of the cbar label)
cbar.set_label(field_label, fontsize=pp.fla, labelpad=pad_cbar_lab)  
    
# Axes labels 
ax.set_xlabel(r'$k_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$y^+$',   fontsize=pp.fla, labelpad=pp.pad_axes_lab)

plt.show()


