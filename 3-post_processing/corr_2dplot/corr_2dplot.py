#!------------------------------------------------------!
#! With this script, we perform plotting of correlation !
#! functions in 2d plots as function of                 !
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

# Import function to setup flow parameters 
from set_flow_parameters import set_flow_parameters

# Import function to plot the correlation functions in 2d plots
from corr_2dplot_sub import corr_2dplot

#!--------------------------------------------------------------------------------------!
#! Main program
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

# Settings for contourf and colormap
cmap_name     = "Greys"
xlabel        = r'$r_z^+$' 
pad_cbar_lab  = -14
size_cbar     = '2%'
    
#!--- Parameters ---!
uwall, nu, twd = set_flow_parameters(itype, re)
 
# Reading of grid points
y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)
y = y[:ny]

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
    
    # Friction Reynolds number
    re_tau = sh_vel * bl_thick / nu
    re_tau = int(re_tau)
    
    # Print friction Reynolds number and boundary layer thickness
    print("Friction Reynolds number, re_tau = ", re_tau)
    print()
    print("Boundary layer thickness, delta_99 = ", bl_thick)
    print()

# Plot 2-dimensional plots for correlation functions
corr_2dplot(Ruuz,'Ruuz',r'$R_{uu}$',Lz,nz,mg_x,nu,y,ny,cmap_name,pad_cbar_lab,size_cbar,add_string,snap_numb,re_tau)
corr_2dplot(Rvvz,'Rvvz',r'$R_{vv}$',Lz,nz,mg_x,nu,y,ny,cmap_name,pad_cbar_lab,size_cbar,add_string,snap_numb,re_tau)
corr_2dplot(Rwwz,'Rwwz',r'$R_{ww}$',Lz,nz,mg_x,nu,y,ny,cmap_name,pad_cbar_lab,size_cbar,add_string,snap_numb,re_tau)
corr_2dplot(Ruvz,'Ruvz',r'$R_{uv}$',Lz,nz,mg_x,nu,y,ny,cmap_name,pad_cbar_lab,size_cbar,add_string,snap_numb,re_tau)







 

    


