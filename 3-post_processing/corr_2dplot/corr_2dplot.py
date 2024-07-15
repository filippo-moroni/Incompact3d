#!------------------------------------------------------!
#! With this script, we perform plotting of correlation !
#! coefficient functions in 2d plots as function of     !
#! y+ and spanwise separation variable rz^+.            !
#!------------------------------------------------------!

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
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
field_label   = r'$R_{uu}^+$'
xlabel        = r'$r_z^+$'
pad_cbar_lab  = -28
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

# Halve the domain size and the number of points in the periodic direction to avoid periodicity effects
# Restriction is imposed also in y if we are dealing with a Channel.
Lz       = Lz / 2.0                       
nz       = nz // 2 
Ruuz     = Ruuz[:ny,:nz]                       
 
# Calculate friction quantities and adimensionalize
sh_vel   = np.sqrt(nu * np.abs(mg_x[0]))  # shear velocity (based on streamwise mean gradient)  
delta_nu = nu / sh_vel                    # viscous length 
y_plus   = y / delta_nu                   
Lz_plus  = Lz / delta_nu
Ruuz     = Ruuz / sh_vel
                                                                                    
# Create the separation variable array
rz = np.linspace(0.0, Lz, nz)
rz = rz / delta_nu

#!--- Plot 1D section ---!

# Limits for axes (used in 'set_plot_settings')
xliminf = 0.0
xlimsup = Lz_plus
yliminf = 0.0
ylimsup = y_plus[-1]

# Extent of the image (dimensions of the domain)
extent = [xliminf, xlimsup, yliminf, ylimsup]

#!--------------------------------------------------------------------------------------!

#!--- Mesh section and iso-levels ---!

# Values of iso-levels        
lvls = np.linspace(np.min(Ruuz), np.max(Ruuz), 20)

field_ticks = [np.min(Ruuz),np.max(Ruuz)]

X, Y = np.meshgrid(rz, y_plus)

#!--------------------------------------------------------------------------------------!

# Auto-correlation function in spanwise direction for longitudinal velocity component

# Subplots environment
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
    
# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Functions to locate the colorbar
divider = make_axes_locatable(ax)
    
cax = divider.append_axes('right', size=size_cbar, pad=0.05)
        
# Imshow function (unexpectedly it adjusts well the aspect ratio of the plotted image with contourf)
im = ax.imshow(Ruuz, cmap=cmap_name, extent=extent, origin='upper')
                
# Plotting with filled contours    
C = ax.contourf(X, Y, Ruuz, lvls, cmap=cmap_name, extend='neither')
    
# Colorbar
cbar = fig.colorbar(C, cax=cax, orientation='vertical', ticks=field_ticks)
       
# Colorbar ticks 
cbar.ax.tick_params(axis='y', labelsize=pp.fla2, length=pp.lmajt, width=pp.tick_width) 
     
# Colorbar label (use pp.pad_cbar_lab to use the default value for padding of the cbar label)
cbar.set_label(field_label, fontsize=pp.fla, labelpad=pad_cbar_lab)  
    
# Axes labels 
ax.set_xlabel(r'$r_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$y^+$',   fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Show the plot
plt.show()


