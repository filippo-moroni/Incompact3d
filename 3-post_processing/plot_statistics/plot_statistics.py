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

# Import functions to read 'input.i3d', 'post.prm' files, statistics data and reference data
from read_files import read_input_files, read_data, read_ref_data

# Import function to read friction Reynolds number Re_tau from .xdmf files
from read_retau import extract_re_tau_value

# Import function to setup flow parameters (kinematic viscosity only at the moment)
from set_flow_parameters import set_flow_parameters

#!--------------------------------------------------------------------------------------!

# Create folders to store later results (e.g. grid spacings and time scales files, plots)
os.makedirs('data_post', mode=0o777, exist_ok=True)
os.makedirs('plots',     mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, numscalar, iswitch_wo, file1, filen, icrfile, nr, add_string = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!
    
#!--- Parameters ---!
uwall, nu, twd = set_flow_parameters(itype, re)

#!--- Reference data ---!
(y_plus_lm,         mean_u_lm, var_u_lm, var_v_lm, mean_uv_lm, 
 rz_plus_cuuz_kim,  cuuz_kim,                          
 rz_plus_cvvz_kim,  cvvz_kim,                          
 rz_plus_cwwz_kim,  cwwz_kim,                          
 y_plus_touber,     mean_u_touber,                     
 y_plus_umean_yao,  mean_u_yao,                        
 y_plus_uvar_yao,   var_u_yao,                         
 y_plus_vvar_yao,   var_v_yao,                         
 y_plus_uvmean_yao, mean_uv_yao,
 y_plus_moser_1999, p_eps_ratio_moser_1999) = read_ref_data() 
 
# Reading of grid points
y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)

# Read statistics data
(mean_u, mean_w, var_u, var_v, mean_uv, 
 vort_x, vort_y, vort_z, mg_tot, mg_x, mg_z,
 eps, Ruuz, Rvvz, Rwwz, Ruvz, Rppz,
 tke_conv, tke_turbt, tke_pstrain, tke_difft, tke_prod, tke_pseps,
 snap_numb) = read_data(itype, numscalar)
                                                                                                                                   
# Valid only for Channel
if itype == 3:
    
    # Change sign for Reynolds stresses
    mean_uv =  - mean_uv

# Valid only for TTBLs
elif itype == 13:
    
    # Initialize the index
    j = 0
    
    # Calculate the index at which the BL thickness delta99 is and delta_99 itself
    while mean_u[j] > mean_u[0]*0.01: 
        
        # Boundary layer thickness delta_99
        bl_thick = y[j]
        
        # Increment the index
        j = j + 1
        
    # Shift due to the translating wall
    mean_u = uwall - mean_u

#!--------------------------------------------------------------------------------------!

#!--- Calculations ---!

# Mesh spacings
delta_x = Lx / nx
delta_z = Lz / nz
           
# Shear quantities
sh_vel   = np.sqrt(nu * np.abs(mg_x[0]))  # shear velocity (based on streamwise mean gradient)  
delta_nu = nu / sh_vel                    # viscous length
t_nu     = nu / (sh_vel ** 2)             # viscous time

# Friction Reynolds number
if itype == 13:
    re_tau = sh_vel * bl_thick / nu
    
    # Print friction Reynolds number and boundary layer thickness
    print("Friction Reynolds number, re_tau = ", re_tau)
    print()
    print("Boundary layer thickness, delta_99 = ", bl_thick)
    print()
    
# Rescaling variables through wall units
delta_x_plus = delta_x / delta_nu
delta_z_plus = delta_z / delta_nu
y_plus       = y       / delta_nu 

Lx_plus = Lx / delta_nu
Ly_plus = Ly / delta_nu 
Lz_plus = Lz / delta_nu

mean_u  /= sh_vel
var_u   /= sh_vel ** 2
var_v   /= sh_vel ** 2
mean_uv /= sh_vel ** 2

# Spanwise velocity is not overwritten since for a channel it is plotted in external units 
mean_w_plus  = mean_w / sh_vel

vort_x *= t_nu
vort_y *= t_nu
vort_z *= t_nu

# Print viscous time unit
print("Viscous time unit, t_nu = ", t_nu)
print()

# y+ at the centerline or at the BL edge
if itype == 3:
    delta_yd_plus = y_plus[ny-1] - y_plus[ny-2]
elif itype == 13:
    delta_yd_plus = y_plus[j] - y_plus[j-1] 

#!--------------------------------------------------------------------------------------!
    
#!--- Writing to file the non-dimensional grid spacings and domain dimensions ---!
          
# Create the file and write  
with open(f'data_post/grid_spacings_post-{snap_numb}_{add_string}.txt', 'w') as f:
    f.write(f"{'delta_x^+':<{pp.c_w}}, "  +
            f"{'delta_y1^+':<{pp.c_w}}, " +
            f"{'delta_z^+':<{pp.c_w}}, "  +
            f"{'Lx^+':<{pp.c_w}}, "       +
            f"{'Ly^+':<{pp.c_w}}, "       +
            f"{'Lz^+':<{pp.c_w}}, "       +
            f"{'delta_yd^+':<{pp.c_w}}\n" )

    f.write(f"{delta_x_plus:{pp.fs}}, "   +
            f"{y_plus[1]:{pp.fs}}, "      +
            f"{delta_z_plus:{pp.fs}}, "   +
            f"{Lx_plus:{pp.fs}}, "        +
            f"{Ly_plus:{pp.fs}}, "        +
            f"{Lz_plus:{pp.fs}}, "        +
            f"{delta_yd_plus:{pp.fs}}\n"  ) 
            
#!--------------------------------------------------------------------------------------!

# Find the maximum of mean total dissipation
eps_max = max(eps)

# Minimum Kolmogorov time scale and print it
tau_eta = np.sqrt(nu/eps_max)
print("Minimum Kolmogorov time scale, tau_eta = ", tau_eta)
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

# Take the correlation functions value at rz = 0 and rescale to obtain correlation coefficients
temp = Ruuz[c,0]
Ruuz = Ruuz / temp

temp = Rvvz[c,0]
Rvvz = Rvvz / temp

temp = Rwwz[c,0]
Rwwz = Rwwz / temp

temp = Ruvz[c,0]
Ruvz = Ruvz / temp

if numscalar == 1:
    temp = Rppz[c,0]
    Rppz = Rppz / temp

# Halve the number of points in z-dir. to avoid periodicity effects
nz = nz // 2
Lz = Lz / 2.0

# Create the separation variable array
rz = np.linspace(0.0, Lz, nz)

# Calculate the index at which the correlation coefficient goes to zero
k = 0
while Ruuz[c,k] > 0.0: k = k + 1

rz0    = rz[0]    # First element of rz vector (rz = 0)
rzstar = rz[k-1]  # Element of rz vector at which Cii(rz) goes to zero

# Rescale separation variable by viscous unit
rz = rz / delta_nu
    
#!--- Calculate the integral length scale lambda z ---!

# Use the minimum between the number of points in the interval Cuu(z) [max, 0]
# and the maximum spline order (5) in order to avoid problem with the spline interpolation itself.
order_spline = min(k-1, 5)

# Interpolation at the 6th order of accuracy with a spline of 5th order
spl = InterpolatedUnivariateSpline(rz[:k], Ruuz[c,:k], k=order_spline)
lambda_z = spl.integral(rz0, rzstar)

# Rescale in wall units
lambda_z = lambda_z / delta_nu

# Print the integral length scale value
print("Integral length scale in viscous units, lambda_z^+ = ", lambda_z)
print()

#!--- Writing to file the viscous time unit and the Kolmogorov time scale ---!
           
# Create the file and write 
with open(f'data_post/time_scales-{snap_numb}_{add_string}.txt', 'w') as f:
    f.write(f"{'t_nu':<{pp.c_w}}, "        +
            f"{'min tau_eta':<{pp.c_w}}\n" )  

    f.write(f"{t_nu:{pp.fs}}, "            +
            f"{tau_eta:{pp.fs}}\n"         )      

#!--------------------------------------------------------------------------------------!

# Ratio between turbulent production and dissipation
p_eps_ratio_tke = np.zeros(ny)

tke_prod  = tke_prod  [:ny]
tke_pseps = tke_pseps [:ny]
p_eps_ratio_tke = np.divide(tke_prod,tke_pseps)

#!--------------------------------------------------------------------------------------!

#!--- Reference mean profiles ---!

# Viscous sub-layer
y_plus_vsl = np.linspace(1, 15, 15)
u_plus_vsl = y_plus_vsl

# Log law constants based on specific flow case
if itype == 13:
   
    # Kozul et al. (2016)
    k = 0.384
    B = 4.173
        
elif itype == 3:

    if pp.iswitch == 0:
    
        # Lee and Moser (2015)
        k = 0.384
        B = 4.27
    
    elif pp.iswitch == 1:
        
        # Cimarelli ('Turbulence' lecture notes)
        k = 0.37
        B = 5.2

# Von Karman law
y_plus_k = np.linspace(5, 180, 175)
u_plus_k = (1.0 / k) * np.log(y_plus_k) + B

#!--------------------------------------------------------------------------------------!
    
#!--- Plot section, mean velocity profile, with selection dipending on the flow case ---!

# Mean velocity profile
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.1
yliminf = 0.0

# Mean velocity profile 
ax.scatter(y_plus[:ny], mean_u[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
# TTBL
if itype == 13:
    
    xlimsup = 520.0
    ylimsup = 30.0
                    
# Channel    
elif itype == 3:

    xlimsup = 300.0
    ylimsup = 25.0
    
    # Lee & Moser (2015)
    ax.plot(y_plus_lm, mean_u_lm, color='C1', linestyle='-', linewidth=pp.lw)
    
    # If wall oscillations are present
    if iswitch_wo == 1:
    
        # Touber & Leschziner (2012)
        #ax.plot(y_plus_touber, mean_u_touber, color='C2', linestyle='-.', linewidth=pp.lw)
        
        # Yao et al. (2019)
        ax.scatter(y_plus_umean_yao, mean_u_yao, marker='^', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='k')
    
# Viscous sublayer and log law
ax.plot(y_plus_vsl, u_plus_vsl, color=pp.grey, linestyle='--', linewidth=pp.lw)
ax.plot(y_plus_k, u_plus_k, color=pp.grey, linestyle='--', linewidth=pp.lw)
    
# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$U^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

# Save and show the figure
save_and_show_plot('umean', snap_numb=snap_numb, add_string=add_string)

#!--------------------------------------------------------------------------------------!

#!--- Mean spanwise velocity profile ---!

# Mean spanwise velocity profile
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# TTBL
if itype == 13:

    # Spanwise mean velocity profile
    ax.scatter(y_plus, mean_w_plus, marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
    # Limits for axes
    xliminf = 0.1
    xlimsup = Ly_plus
    yliminf = min(mean_w_plus)*1.2    
    ylimsup = max(mean_w_plus)*1.2
    
    # Axes labels
    ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$W^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

# Channel
elif itype == 3:

    # Spanwise mean velocity profile
    ax.scatter(y[:ny], mean_w[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
    # Limits for axes
    xliminf = 0.0
    xlimsup = 1.0
    yliminf = min(mean_w)*1.2    
    ylimsup = max(mean_w)*1.2
    
    # Axes labels
    ax.set_xlabel(r'$y/h$',   fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$W/U_p$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Save and show the figure
save_and_show_plot('wmean', snap_numb=snap_numb, add_string=add_string)

#!--------------------------------------------------------------------------------------!

# <u'u'>
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.1
yliminf = 0.0

# <u'u'>
ax.scatter(y_plus[:ny], var_u[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

# TTBL
if itype == 13:
    
    xlimsup = 520.0
    ylimsup = 10.0       
           
# Channel    
elif itype == 3:

    xlimsup = 300.0
    ylimsup = 8.0
    
    # Lee & Moser (2015)
    ax.plot(y_plus_lm, var_u_lm, color='C1', linestyle='-', linewidth=pp.lw)
    
    # If wall oscillations are present
    if iswitch_wo == 1:
    
        # Yao et al. (2019)
        ax.scatter(y_plus_uvar_yao, var_u_yao, marker='^', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='k')
    
# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$\langle u^{\prime 2} \rangle^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

# Save and show the figure
save_and_show_plot('uvar', snap_numb=snap_numb, add_string=add_string)

#!--------------------------------------------------------------------------------------!

# <v'v'>
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.1
yliminf = 0.0

# <v'v'>
ax.scatter(y_plus[:ny], var_v[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

# TTBL
if itype == 13:

    xlimsup = 520.0
    ylimsup = 1.2
            
# Channel    
elif itype == 3:

    xlimsup = 300.0
    ylimsup = 0.8
    
    # Lee & Moser (2015)
    ax.plot(y_plus_lm, var_v_lm, color='C1', linestyle='-', linewidth=pp.lw)
    
    # If wall oscillations are present
    if iswitch_wo == 1:
    
        # Yao et al. (2019)
        ax.scatter(y_plus_vvar_yao, var_v_yao, marker='^', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='k')
    
# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$\langle v^{\prime 2} \rangle^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

# Save and show the figure
save_and_show_plot('vvar', snap_numb=snap_numb, add_string=add_string)

#!--------------------------------------------------------------------------------------!

# <u'v'>
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.1
yliminf = 0.0

# <u'v'>
ax.scatter(y_plus[:ny], mean_uv[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
# TTBL
if itype == 13:

    xlimsup = 520.0
    ylimsup = 1.2
        
    # y-axis label
    ax.set_ylabel(r'$\langle u^{\prime} v^{\prime}\rangle^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
        
# Channel    
elif itype == 3:

    xlimsup = 300.0
    ylimsup = 0.8
        
    # Lee & Moser (2015)
    ax.plot(y_plus_lm, mean_uv_lm, color='C1', linestyle='-', linewidth=pp.lw)
       
    # y-axis label
    ax.set_ylabel(r'$-\langle u^{\prime} v^{\prime}\rangle^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    # If wall oscillations are present
    if iswitch_wo == 1:
    
        # Yao et al. (2019)
        ax.scatter(y_plus_uvmean_yao, mean_uv_yao, marker='^', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='k')
 
# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

# Save and show the figure
save_and_show_plot('uvmean', snap_numb=snap_numb, add_string=add_string)

#!--------------------------------------------------------------------------------------!

# Correlation coefficients in spanwise direction (z)

#!--------------------------------------------------------------------------------------!

# Cuuz
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.0
xlimsup = Lz_plus / 2.0
yliminf = np.min(Ruuz[c,:])*1.2
ylimsup = np.max(Ruuz[c,:])*1.2

# Auto-correlation coefficient for u'
ax.scatter(rz, Ruuz[c,:nz], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

# Reference data at y+ = 10 
if y_plus_in == 10.0:

    # Kim et al. (1987) data
    ax.plot(rz_plus_cuuz_kim, cuuz_kim, color='C1', linestyle='-', linewidth=pp.lw)

# Plot horizontal line at Cuu = 0
ax.hlines(y=0.0, xmin=xliminf, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed')
                    
# Axes labels
ax.set_xlabel(r'$r_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$C_{uu}(r_z^+)$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Save and show the figure
save_and_show_plot('Cuuz', snap_numb=snap_numb, add_string=add_string, y_plus_in=y_plus_in)

#!--------------------------------------------------------------------------------------!

# Cvvz
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.0
xlimsup = Lz_plus / 2.0
yliminf = np.min(Rvvz[c,:])*1.2
ylimsup = np.max(Rvvz[c,:])*1.2

# Auto-correlation coefficient for v'
ax.scatter(rz, Rvvz[c,:nz], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

# Reference data at y+ = 10 
if y_plus_in == 10.0:

    # Kim et al. (1987) data
    ax.plot(rz_plus_cvvz_kim, cvvz_kim, color='C1', linestyle='-', linewidth=pp.lw)

# Plot horizontal line at Cvv = 0
ax.hlines(y=0.0, xmin=xliminf, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed')
                  
# Axes labels
ax.set_xlabel(r'$r_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$C_{vv}(r_z^+)$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Save and show the figure
save_and_show_plot('Cvvz', snap_numb=snap_numb, add_string=add_string, y_plus_in=y_plus_in)

#!--------------------------------------------------------------------------------------!

# Cwwz
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.0
xlimsup = Lz_plus / 2.0
yliminf = np.min(Rwwz[c,:])*1.2
ylimsup = np.max(Rwwz[c,:])*1.2

# Auto-correlation coefficient for w'
ax.scatter(rz, Rwwz[c,:nz], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

# Reference data at y+ = 10 
if y_plus_in == 10.0:

    # Kim et al. (1987) data
    ax.plot(rz_plus_cwwz_kim, cwwz_kim, color='C1', linestyle='-', linewidth=pp.lw)

# Plot horizontal line at Cww = 0
ax.hlines(y=0.0, xmin=xliminf, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed')
                   
# Axes labels
ax.set_xlabel(r'$r_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$C_{ww}(r_z^+)$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Save and show the figure
save_and_show_plot('Cwwz', snap_numb=snap_numb, add_string=add_string, y_plus_in=y_plus_in)

#!--------------------------------------------------------------------------------------!

# Cuvz
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.0
xlimsup = Lz_plus / 2.0
yliminf = np.min(Ruvz[c,:])*1.2
ylimsup = np.max(Ruvz[c,:])*1.2

# Correlation coefficient for u' and v'
ax.scatter(rz, Ruvz[c,:nz], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

# Plot horizontal line at Cuv = 0
ax.hlines(y=0.0, xmin=xliminf, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed')
                   
# Axes labels
ax.set_xlabel(r'$r_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$C_{uv}(r_z^+)$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Save and show the figure
save_and_show_plot('Cuvz', snap_numb=snap_numb, add_string=add_string, y_plus_in=y_plus_in)

#!--------------------------------------------------------------------------------------!

if numscalar == 1:

    # Cppz
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

    # Limits for axes
    xliminf = 0.0
    xlimsup = Lz_plus / 2.0
    yliminf = np.min(Rppz[c,:])*1.2
    ylimsup = np.max(Rppz[c,:])*1.2

    # Auto-correlation coefficient for phi'
    ax.scatter(rz, Rppz[c,:nz], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

    # Plot horizontal line at Cpp = 0
    ax.hlines(y=0.0, xmin=xliminf, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed')
                   
    # Axes labels
    ax.set_xlabel(r'$r_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$C_{\varphi \varphi}(r_z^+)$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

    # Save and show the figure
    save_and_show_plot('Cppz', snap_numb=snap_numb, add_string=add_string, y_plus_in=y_plus_in)

#!--------------------------------------------------------------------------------------!

# Ratio of production over dissipation of TKE
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.1
yliminf = 0.0
ylimsup = 2.0

# P/eps of TKE
ax.scatter(y_plus[:ny], p_eps_ratio_tke[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
# TTBL
if itype == 13:

    xlimsup = 520.0
                
# Channel    
elif itype == 3:

    xlimsup = 300.0
        
    # Moser et al. (1999)
    ax.plot(y_plus_moser_1999, p_eps_ratio_moser_1999, color='C1', linestyle='-', linewidth=pp.lw)
            
# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$P/\varepsilon$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Save and show the figure
save_and_show_plot('p_eps_ratio_tke', snap_numb=snap_numb, add_string=add_string)

#!--------------------------------------------------------------------------------------!



