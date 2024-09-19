
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this script, we perform plotting of spanwise correlation 
!              functions in 2d plots as function of y+ and 
!              spanwise separation variable rz^+.            
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
import sys
import os
import numpy as np

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '..', 'python_common'))
sys.path.append(config_path)

# Import functions to read 'input.i3d', 'post.prm' files and statistics data
from read_files import read_input_files, read_data

# Import function to setup flow parameters 
from set_flow_parameters import set_flow_parameters

# Import function to plot the correlation functions in 2d plots
from corr_2dplot_sub import corr_2dplot

# Import function to calculate boundary layer thickness delta_99 for a TTBL
from ttbl_subs import calculate_ttbl_delta_99

#!--------------------------------------------------------------------------------------!
#! Main program
#!--------------------------------------------------------------------------------------!

# Create folders to store later results (e.g. grid spacings and time scales files, plots)
os.makedirs('plots',     mode=0o777, exist_ok=True)
os.makedirs('plots/correlations',        mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
(itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, ioutput, iswitch_wo,  
 add_string, file1, filen, icrfile, nr, post_mean, post_vort, post_diss, post_corz, post_tke_eq
) = read_input_files('input.i3d','post.prm')

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
 eps, Ruuz, Rvvz, Rwwz, Ruvz, Rssz,
 tke_turbt, tke_presst, tke_difft, tke_prod, tke_pseps,
 snap_numb) = read_data(itype, numscalar, post_mean, post_vort, post_diss, 
                        post_corz, post_tke_eq, ny, nz)

# Inner quantities
sh_vel   = np.sqrt(nu * np.abs(mg_x[0]))  # shear velocity (based on streamwise mean gradient)  

# Valid only for TTBLs
if itype == 13:
    
    # Calculate BL thickness delta_99 for a TTBL and its related index
    (bl_thick, bl_thick_j) = calculate_ttbl_delta_99(mean_u, y)
    
    # Friction Reynolds number
    re_tau = sh_vel * bl_thick / nu
    re_tau = int(re_tau)
    
    # Print friction Reynolds number and boundary layer thickness
    print(">>> Friction Reynolds number, re_tau = ", re_tau)
    print()
    print(">>> Boundary layer thickness, delta_99 = ", bl_thick)
    print()

# Plot 2-dimensional plots for correlation functions
corr_2dplot(Ruuz,'Ruuz',r'$R_{uu}$',Lz,nz,mg_x,nu,y,ny,cmap_name,pad_cbar_lab,size_cbar,add_string,snap_numb,re_tau)
corr_2dplot(Rvvz,'Rvvz',r'$R_{vv}$',Lz,nz,mg_x,nu,y,ny,cmap_name,pad_cbar_lab,size_cbar,add_string,snap_numb,re_tau)
corr_2dplot(Rwwz,'Rwwz',r'$R_{ww}$',Lz,nz,mg_x,nu,y,ny,cmap_name,pad_cbar_lab,size_cbar,add_string,snap_numb,re_tau)
corr_2dplot(Ruvz,'Ruvz',r'$R_{uv}$',Lz,nz,mg_x,nu,y,ny,cmap_name,pad_cbar_lab,size_cbar,add_string,snap_numb,re_tau)

if numscalar == 1:
    corr_2dplot(Rssz,'Rssz',r'$R_{\varphi\varphi}$',Lz,nz,mg_x,nu,y,ny,cmap_name,pad_cbar_lab,size_cbar,add_string,snap_numb,re_tau)





 

    


