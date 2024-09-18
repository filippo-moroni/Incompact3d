
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this script, we perform plotting of statistics for TTBLs 
!              and channel flow simulations:                 
!               
!               - mean statistics (mean[u], mean[w], 
!                                  var[u], var[v], var[w],
!                                  mean[uv]);              
!               - mean total dissipation;                  
!               - correlation coefficients (Cuu, Cvv, Cww, Cuv, Css) for 
!                 spanwise correlations;
!               - Turbulent Kinetic Energy (TKE) budget and production
!                 over dissipation ratio.  
!                                                         
!              Calculated and stored:                                  
!                                                         
!               - non-dimensional grid spacings and domain dimensions;
!               - time-scales (viscous time and Kolmogorov time scale).
!                         
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
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

# Import functions to setting up, save and show plots and to get 
# reference mean streamwise velocity profile  
from plot_subs import set_plot_settings, save_and_show_plot, get_ref_mean_vel_profile

# Import functions to read 'input.i3d', 'post.prm' files, statistics data and reference data
from read_files import read_input_files, read_data, read_ref_data

# Import function to setup flow parameters 
from set_flow_parameters import set_flow_parameters

#!--------------------------------------------------------------------------------------!

# Create folders to store later results (e.g. grid spacings and time scales files, plots)
os.makedirs('data_post',                 mode=0o777, exist_ok=True)
os.makedirs('data_post/plot_statistics', mode=0o777, exist_ok=True)
os.makedirs('plots',                     mode=0o777, exist_ok=True)
os.makedirs('plots/mean_stats',          mode=0o777, exist_ok=True)
os.makedirs('plots/vort_stats',          mode=0o777, exist_ok=True)
os.makedirs('plots/diss_stats',          mode=0o777, exist_ok=True)
os.makedirs('plots/correlations',        mode=0o777, exist_ok=True)
os.makedirs('plots/tke_stats',           mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
(itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, ioutput, iswitch_wo,  
 add_string, file1, filen, icrfile, nr, post_mean, post_vort, post_diss, post_corz, post_tke_eq
) = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!
    
#!--- Parameters ---!
uwall, nu, twd = set_flow_parameters(itype, re)

#!--- Reference data ---!
(y_plus_lm,                 mean_u_lm, var_u_lm, var_v_lm, var_w_lm, mean_uv_lm,
 rz_plus_cuuz_kim,          cuuz_kim, 
 rz_plus_cvvz_kim,          cvvz_kim,
 rz_plus_cwwz_kim,          cwwz_kim,
 y_plus_touber,             mean_u_touber,
 y_plus_umean_yao,          mean_u_yao,
 y_plus_uvar_yao,           var_u_yao,
 y_plus_vvar_yao,           var_v_yao,
 y_plus_wvar_yao,           var_w_yao,
 y_plus_uvmean_yao,         mean_uv_yao, 
 y_plus_moser_1999,         p_eps_ratio_moser_1999,
 y_plus_lm1000,             p_eps_ratio_lm1000,
 y_plus_umean_kozul,        mean_u_kozul,
 y_plus_uvar_kozul,         var_u_kozul,
 y_plus_vvar_kozul,         var_v_kozul,
 y_plus_uvmean_kozul,       mean_uv_kozul,
 y_plus_tke_turbt_mansour,  tke_turbt_mansour,  
 y_plus_tke_presst_mansour, tke_presst_mansour,
 y_plus_tke_difft_mansour,  tke_difft_mansour,
 y_plus_tke_prod_mansour,   tke_prod_mansour,
 y_plus_tke_pseps_mansour,  tke_pseps_mansour) = read_ref_data()
  
# Reading of grid points
y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)

# Read statistics data
(mean_u, mean_w, var_u, var_v, var_w, mean_uv, 
 vort_x, vort_y, vort_z, mg_tot, mg_x, mg_z,
 eps, Ruuz, Rvvz, Rwwz, Ruvz, Rssz,
 tke_turbt, tke_presst, tke_difft, tke_prod, tke_pseps,
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
  
# Rescaling variables through wall units
delta_x_plus = delta_x / delta_nu
delta_z_plus = delta_z / delta_nu
y_plus       = y       / delta_nu 

Lx_plus = Lx / delta_nu
Ly_plus = Ly / delta_nu 
Lz_plus = Lz / delta_nu

# Rescale mean flow statistics
if post_mean:
    mean_u  /= sh_vel
    var_u   /= sh_vel**2
    var_v   /= sh_vel**2
    var_w   /= sh_vel**2
    mean_uv /= sh_vel**2

    # Spanwise velocity is not overwritten since for a channel it is plotted in external units 
    mean_w_plus  = mean_w / sh_vel

# Rescale vorticity
if post_vort:
    vort_x *= t_nu
    vort_y *= t_nu
    vort_z *= t_nu

# Rescale TKE terms
if post_tke_eq:
    tke_turbt  /= sh_vel**2
    tke_presst /= sh_vel**2
    tke_difft  /= sh_vel**2
    tke_prod   /= sh_vel**2
    tke_pseps  /= sh_vel**2

# Friction Reynolds number
re_tau = None

if itype == 13:
    re_tau = sh_vel * bl_thick / nu
    
    # Remove decimal digits
    re_tau = int(re_tau)
    
    # Print friction Reynolds number and boundary layer thickness
    print(">>> Friction Reynolds number, Re_tau = ", re_tau)
    print()
    print(">>> Boundary layer thickness, delta_99 = ", bl_thick)
    print()
    print(">>> Domain height in wall units, Ly+ = ", Ly_plus)
    print()

# Print viscous time unit
print(">>> Viscous time unit, t_nu = ", t_nu)
print()

# y+ at the centerline or at the BL edge and halving Ly+ in case of Channel
if itype == 3:
    delta_yd_plus = y_plus[ny-1] - y_plus[ny-2]
    Ly_plus = Ly_plus / 2.0
elif itype == 13:
    delta_yd_plus = y_plus[j] - y_plus[j-1] 

#!--------------------------------------------------------------------------------------!
    
#!--- Writing to file the non-dimensional grid spacings and domain dimensions ---!

print(">>> Saving in 'grid_spacings_post' non-dimensional grid spacings")
print(">>> and domain dimensions.")
print(">>> Folder: data_post/plot_statistics/.")
print()

# Create the file and write  
with open(f'data_post/plot_statistics/grid_spacings_post-{snap_numb}_{add_string}.txt', 'w') as f:
    f.write(f"{'delta_x^+':>{pp.c_w}}, "  +
            f"{'delta_y1^+':>{pp.c_w}}, " +
            f"{'delta_z^+':>{pp.c_w}}, "  +
            f"{'Lx^+':>{pp.c_w}}, "       +
            f"{'Ly^+':>{pp.c_w}}, "       +
            f"{'Lz^+':>{pp.c_w}}, "       +
            f"{'delta_yd^+':>{pp.c_w}}\n" )

    f.write(f"{delta_x_plus:{pp.fs}}, "   +
            f"{y_plus[1]:{pp.fs}}, "      +
            f"{delta_z_plus:{pp.fs}}, "   +
            f"{Lx_plus:{pp.fs}}, "        +
            f"{Ly_plus:{pp.fs}}, "        +
            f"{Lz_plus:{pp.fs}}, "        +
            f"{delta_yd_plus:{pp.fs}}\n"  ) 
            
#!--------------------------------------------------------------------------------------!

#!--- Plotting mean statistics ---!
if post_mean:

    print(">>> Plotting mean streamwise velocity profile.")
    print(">>> Folder: plots/mean_stats/.")
    print()

    # Get reference mean streamwise velocity profile
    (y_plus_vsl,u_plus_vsl,y_plus_k,u_plus_k) = get_ref_mean_vel_profile(itype,pp.iswitch)

    # Mean velocity profile
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

    # Limits for axes
    xliminf = 0.1
    yliminf = 0.0

    # Mean velocity profile 
    ax.scatter(y_plus[:ny], mean_u[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
    # Viscous sublayer and log law
    ax.plot(y_plus_vsl, u_plus_vsl, color=pp.grey, linestyle='--', linewidth=pp.lw)
    ax.plot(y_plus_k,   u_plus_k,   color=pp.grey, linestyle='--', linewidth=pp.lw)
    
    # Axes labels
    ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$U^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    # Description of .pdf file
    description = 'Mean streamwise velocity profile.'
    
    # TTBL
    if itype == 13:
    
        xlimsup = Ly_plus
        ylimsup = 30.0
     
        # Plot Lee & Moser (2015) data if Re_tau is close to their Re_tau (180) 
        if 150 < re_tau < 210:
        
            # Lee & Moser (2015)
            ax.plot(y_plus_lm, mean_u_lm, color='C1', linestyle='-', linewidth=pp.lw)
            
            # Completing description
            description += ' Reference data Lee & Moser (2015).'
                        
        # Plot Kozul et al. (2016) data if Re_tau is close to their Re_tau (~ 425) (Re_theta = 1100)
        elif 400 < re_tau < 450:
        
            # Kozul et al. (2016)
            ax.plot(y_plus_umean_kozul, mean_u_kozul, color='C1', linestyle='-', linewidth=pp.lw)
            
            # Completing description
            description += ' Reference data Kozul et al. (2016).'
                           
    # Channel    
    elif itype == 3:

        xlimsup = 300.0
        ylimsup = 25.0
    
        # Lee & Moser (2015)
        ax.plot(y_plus_lm, mean_u_lm, color='C1', linestyle='-', linewidth=pp.lw)
        
        # Completing description
        description += ' Reference data Lee & Moser (2015).'
    
        # If wall oscillations are present
        if iswitch_wo == 1:
    
            # Touber & Leschziner (2012)
            #ax.plot(y_plus_touber, mean_u_touber, color='C2', linestyle='-.', linewidth=pp.lw)
        
            # Yao et al. (2019)
            ax.scatter(y_plus_umean_yao, mean_u_yao, marker='^', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='k')
            
            # Completing description, with slicing of the string (removing the final dot)
            description = description[:-1]
            
            description += ', Yao et al. (2019).'
    
    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)
    
    # Save and show the figure
    save_and_show_plot('umean', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, subfolder='mean_stats', description=description)

    #!--------------------------------------------------------------------------------------!

    print(">>> Plotting mean spanwise velocity profile.")
    print(">>> Folder: plots/mean_stats/.")
    print()

    # Mean spanwise velocity profile
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
    
    # Description of .pdf file
    description = 'Mean spanwise velocity profile.'

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
        ax.scatter(y, mean_w, marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
        # Limits for axes
        xliminf = 0.0
        xlimsup = 2.0
        yliminf = min(mean_w)*1.2    
        ylimsup = max(mean_w)*1.2
    
        # Axes labels
        ax.set_xlabel(r'$y/h$',   fontsize=pp.fla, labelpad=pp.pad_axes_lab)
        ax.set_ylabel(r'$W/U_p$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
        # Set the plot parameters using the function 'set_plot_settings'
        # Last argument is the switcher for semilog plot (1: yes, 0: no)
        set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

    # Save and show the figure
    save_and_show_plot('wmean', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, subfolder='mean_stats', description=description)
    
    #!--------------------------------------------------------------------------------------!

    print(">>> Plotting streamwise velocity variance.")
    print(">>> Folder: plots/mean_stats/.")
    print()

    # <u'u'>
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

    # Limits for axes
    xliminf = 0.1
    yliminf = 0.0

    # <u'u'>
    ax.scatter(y_plus[:ny], var_u[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
    # Description of .pdf file
    description = 'Streamwise velocity variance.'

    # TTBL
    if itype == 13:
    
        xlimsup = Ly_plus
        ylimsup = 10.0
        
        # Plot Lee & Moser (2015) data if Re_tau is close to their Re_tau (180) 
        if 150 < re_tau < 210:
        
            # Lee & Moser (2015)
            ax.plot(y_plus_lm, var_u_lm, color='C1', linestyle='-', linewidth=pp.lw)
            
            # Completing description
            description += ' Reference data Lee & Moser (2015).'
        
        # Plot Kozul et al. (2016) data if Re_tau is close to their Re_tau (~ 425) (Re_theta = 1100)
        elif 400 < re_tau < 450:
        
            # Kozul et al. (2016)
            ax.plot(y_plus_uvar_kozul, var_u_kozul, color='C1', linestyle='-', linewidth=pp.lw)
            
            # Completing description
            description += ' Reference data Kozul et al. (2016).'    
               
    # Channel    
    elif itype == 3:

        xlimsup = 300.0
        ylimsup = 8.0
    
        # Lee & Moser (2015)
        ax.plot(y_plus_lm, var_u_lm, color='C1', linestyle='-', linewidth=pp.lw)
        
        # Completing description
        description += ' Reference data Lee & Moser (2015).'
    
        # If wall oscillations are present
        if iswitch_wo == 1:
    
            # Yao et al. (2019)
            ax.scatter(y_plus_uvar_yao, var_u_yao, marker='^', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='k')
            
            # Completing description, with slicing of the string (removing the final dot)
            description = description[:-1]
            
            description += ', Yao et al. (2019).'
    
    # Axes labels
    ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$\langle u^{\prime 2} \rangle^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

    # Save and show the figure
    save_and_show_plot('uvar', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, subfolder='mean_stats', description=description)

    #!--------------------------------------------------------------------------------------!

    print(">>> Plotting wall-normal velocity variance.")
    print(">>> Folder: plots/mean_stats/.")
    print()

    # <v'v'>
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

    # Limits for axes
    xliminf = 0.1
    yliminf = 0.0

    # <v'v'>
    ax.scatter(y_plus[:ny], var_v[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
    # Description of .pdf file
    description = 'Wall-normal velocity variance.'

    # TTBL
    if itype == 13:

        xlimsup = Ly_plus
        ylimsup = 1.5
        
        # Plot Lee & Moser (2015) data if Re_tau is close to their Re_tau (180) 
        if 150 < re_tau < 210:
        
            # Lee & Moser (2015)
            ax.plot(y_plus_lm, var_v_lm, color='C1', linestyle='-', linewidth=pp.lw)
            
            # Completing description
            description += ' Reference data Lee & Moser (2015).'
            
        # Plot Kozul et al. (2016) data if Re_tau is close to their Re_tau (~ 425) (Re_theta = 1100)
        elif 400 < re_tau < 450:
        
            # Kozul et al. (2016)
            ax.plot(y_plus_vvar_kozul, var_v_kozul, color='C1', linestyle='-', linewidth=pp.lw)
            
            # Completing description
            description += ' Reference data Kozul et al. (2016).'  
            
    # Channel    
    elif itype == 3:

        xlimsup = 300.0
        ylimsup = 0.8
    
        # Lee & Moser (2015)
        ax.plot(y_plus_lm, var_v_lm, color='C1', linestyle='-', linewidth=pp.lw)
        
        # Completing description
        description += ' Reference data Lee & Moser (2015).'
    
        # If wall oscillations are present
        if iswitch_wo == 1:
    
            # Yao et al. (2019)
            ax.scatter(y_plus_vvar_yao, var_v_yao, marker='^', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='k')
            
            # Completing description, with slicing of the string (removing the final dot)
            description = description[:-1]
            
            description += ', Yao et al. (2019).'
    
    # Axes labels
    ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$\langle v^{\prime 2} \rangle^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

    # Save and show the figure
    save_and_show_plot('vvar', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, subfolder='mean_stats', description=description)

    #!--------------------------------------------------------------------------------------!
    
    print(">>> Plotting spanwise velocity variance.")
    print(">>> Folder: plots/mean_stats/.")
    print()

    # <w'w'>
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

    # Limits for axes
    xliminf = 0.1
    yliminf = 0.0

    # <w'w'>
    ax.scatter(y_plus[:ny], var_w[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
    # Description of .pdf file
    description = 'Spanwise velocity variance.'

    # TTBL
    if itype == 13:

        xlimsup = Ly_plus
        ylimsup = 2.5
        
        # Plot Lee & Moser (2015) data if Re_tau is close to their Re_tau (180) 
        if 150 < re_tau < 210:
        
            # Lee & Moser (2015)
            ax.plot(y_plus_lm, var_w_lm, color='C1', linestyle='-', linewidth=pp.lw)
            
            # Completing description
            description += ' Reference data Lee & Moser (2015).'
            
    # Channel    
    elif itype == 3:

        xlimsup = 300.0
        ylimsup = 1.5
    
        # Lee & Moser (2015)
        ax.plot(y_plus_lm, var_w_lm, color='C1', linestyle='-', linewidth=pp.lw)
        
        # Completing description
        description += ' Reference data Lee & Moser (2015).'
    
        # If wall oscillations are present
        if iswitch_wo == 1:
    
            # Yao et al. (2019)
            ax.scatter(y_plus_wvar_yao, var_w_yao, marker='^', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='k')
            
            # Completing description, with slicing of the string (removing the final dot)
            description = description[:-1]
            
            description += ', Yao et al. (2019).'
    
    # Axes labels
    ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$\langle w^{\prime 2} \rangle^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

    # Save and show the figure
    save_and_show_plot('wvar', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, subfolder='mean_stats', description=description)

    #!--------------------------------------------------------------------------------------!

    print(">>> Plotting Reynolds stress <u'v'>.")
    print(">>> Folder: plots/mean_stats/.")
    print()

    # <u'v'>
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

    # Limits for axes
    xliminf = 0.1
    yliminf = 0.0

    # <u'v'>
    ax.scatter(y_plus[:ny], mean_uv[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
    # Description of .pdf file
    description = 'Reynolds stress.'
    
    # TTBL
    if itype == 13:

        xlimsup = Ly_plus
        ylimsup = 1.2
        
        # y-axis label
        ax.set_ylabel(r'$\langle u^{\prime} v^{\prime}\rangle^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
        
        # Plot Lee & Moser (2015) data if Re_tau is close to their Re_tau (180) 
        if 150 < re_tau < 210:
        
            # Lee & Moser (2015)
            ax.plot(y_plus_lm, mean_uv_lm, color='C1', linestyle='-', linewidth=pp.lw)
            
            # Completing description
            description += ' Reference data Lee & Moser (2015).'
            
        # Plot Kozul et al. (2016) data if Re_tau is close to their Re_tau (~ 425) (Re_theta = 1100)
        elif 400 < re_tau < 450:
        
            # Kozul et al. (2016)
            ax.plot(y_plus_uvmean_kozul, mean_uv_kozul, color='C1', linestyle='-', linewidth=pp.lw)
            
            # Completing description
            description += ' Reference data Kozul et al. (2016).'   
            
    # Channel    
    elif itype == 3:

        xlimsup = 300.0
        ylimsup = 0.8
        
        # Lee & Moser (2015)
        ax.plot(y_plus_lm, mean_uv_lm, color='C1', linestyle='-', linewidth=pp.lw)
        
        # Completing description
        description += ' Reference data Lee & Moser (2015).'
       
        # y-axis label
        ax.set_ylabel(r'$-\langle u^{\prime} v^{\prime}\rangle^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
        # If wall oscillations are present
        if iswitch_wo == 1:
    
            # Yao et al. (2019)
            ax.scatter(y_plus_uvmean_yao, mean_uv_yao, marker='^', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='k')
            
            # Completing description, with slicing of the string (removing the final dot)
            description = description[:-1]
            
            description += ', Yao et al. (2019).'
 
    # Axes labels
    ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

    # Save and show the figure
    save_and_show_plot('uvmean', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, subfolder='mean_stats', description=description)

    #!--------------------------------------------------------------------------------------!

#!--- Total dissipation section ---!
if post_diss:

    # Find the maximum of mean total dissipation
    eps_max = max(eps)

    # Minimum Kolmogorov time scale and print it
    tau_eta = np.sqrt(nu/eps_max)
    print("Minimum Kolmogorov time scale, tau_eta = ", tau_eta)
    print()

    #!--- Writing to file the viscous time unit and the Kolmogorov time scale ---!
    print(">>> Saving in 'time_scales' viscous time unit and")
    print(">>> Kolmogorov time scale.")
    print(">>> Folder: data_post/plot_statistics/.")
    print()
            
    # Create the file and write 
    with open(f'data_post/plot_statistics/time_scales-{snap_numb}_{add_string}.txt', 'w') as f:
        f.write(f"{'t_nu':>{pp.c_w}}, "        +
                f"{'min tau_eta':>{pp.c_w}}\n" )  

        f.write(f"{t_nu:{pp.fs}}, "            +
                f"{tau_eta:{pp.fs}}\n"         )

    print(">>> Plotting total dissipation.")
    print(">>> Folder: plots/diss_stats/.")
    print()

    # Total dissipation
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
    
    # Description of .pdf file
    description = 'Total dissipation.'

    # Total dissipation
    ax.scatter(y_plus, eps, marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
    # Limits for axes
    xliminf = 0.1
    xlimsup = Ly_plus
    yliminf = min(eps)*1.2    
    ylimsup = max(eps)*1.2
    
    # Axes labels
    ax.set_xlabel(r'$y^+$',               fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$\varepsilon_{tot}$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

    # Save and show the figure
    save_and_show_plot('eps_tot', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, subfolder='diss_stats', description=description)

    #!--------------------------------------------------------------------------------------!

#!--- Correlation coefficients in spanwise direction (z) ---!
if post_corz:

    #!--- Calculations for correlations ---!

    # Select the height at which correlations are plotted
    y_plus_in = np.float64(input(">>> Enter y+ value for correlations plotting: "))
    print()

    # Search for the index corresponding to the target y+ for correlations
    c = 0 
    for j in range(0, ny-1, 1):   
        if y_plus[j] < y_plus_in: c = c + 1

    # Print the actual y+ value selected
    print(">>> Actual y+ value selected = ", y_plus[c])
    print()

    # Print the corresponding j-th index
    print(">>> Corresponding j-th index = ", c)
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
        temp = Rssz[c,0]
        Rssz = Rssz / temp

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
    print(">>> Integral length scale in viscous units, lambda_z^+ = ", lambda_z)
    print()
    
    #!--------------------------------------------------------------------------------------!

    # Define data and labels for each auto-correlation coefficient Cii(rz^+) plot
    cii_plots = [
                 {
                  'title': 'Cuuz',
                  'var': Ruuz,
                  'kim_data': cuuz_kim,
                  'kim_rz_plus': rz_plus_cuuz_kim,
                  'ylabel': r'$C_{uu}(r_z^+)$',
                  
                  # 'Local' description
                  'descr': 'Auto-correlation coefficient for streamwise velocity component in spanwise direction.' 
                 },
                 {
                  'title': 'Cvvz',
                  'var': Rvvz,
                  'kim_data': cvvz_kim,
                  'kim_rz_plus': rz_plus_cvvz_kim,
                  'ylabel': r'$C_{vv}(r_z^+)$',
                  
                  # 'Local' description
                  'descr': 'Auto-correlation coefficient for wall-normal velocity component in spanwise direction.' 
                 },
                 {
                  'title': 'Cwwz',
                  'var': Rwwz,
                  'kim_data': cwwz_kim,
                  'kim_rz_plus': rz_plus_cwwz_kim,
                  'ylabel': r'$C_{ww}(r_z^+)$',
                  
                  # 'Local' description
                  'descr': 'Auto-correlation coefficient for spanwise velocity component in spanwise direction.' 
                 }
                ]
    
    print(">>> Plotting spanwise velocity auto-correlation coefficients Cii(rz).")
    print(">>> Folder: plots/correlations/.")
    print()

    # Loop through each plot's data and create the plots
    for plot in cii_plots:
        
        # Create the figure and axis
        fig, ax = plt.subplots(1, 1, figsize=(pp.xinches, pp.yinches), linewidth=pp.tick_width, dpi=300)

        # Limits for axes
        xliminf = 0.0
        xlimsup = Lz_plus / 2.0

        min_value1 = np.min(plot['var'][c, :])
        min_value2 = np.min(plot['kim_data'][:])

        yliminf = min(min_value1, min_value2) * 1.2
        ylimsup = 1.2

        # Auto-correlation coefficient Cii(rz^+)
        ax.scatter(rz, plot['var'][c, :nz], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

        # Reference data at y+ = 10
        if y_plus_in == 10.0:
    
            # Kim et al. (1987) data
            ax.plot(plot['kim_rz_plus'], plot['kim_data'], color='C1', linestyle='-', linewidth=pp.lw)
            
            # Complete description            
            plot['descr'] += ' Reference data Kim et al. (1987).'
            

        # Plot horizontal line at C = 0
        ax.hlines(y=0.0, xmin=xliminf, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed')

        # Axes labels
        ax.set_xlabel(r'$r_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
        ax.set_ylabel(plot['ylabel'], fontsize=pp.fla, labelpad=pp.pad_axes_lab)

        # Set the plot parameters using the function 'set_plot_settings'
        # Last argument is the switcher for semilog plot (1: yes, 0: no)
        set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)
        
        # Save and show the figure
        save_and_show_plot(plot['title'], snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, y_plus_in=y_plus_in, subfolder='correlations', description=plot['descr'])

    #!--------------------------------------------------------------------------------------!

    print(">>> Plotting spanwise velocity correlation coefficient Cuv(rz).")
    print(">>> Folder: plots/correlations/.")
    print()

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
    
    # Description of .pdf file
    description = 'Correlation coefficient for streamwise and wall-normal velocity components in spanwise direction.'

    # Save and show the figure
    save_and_show_plot('Cuvz', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, y_plus_in=y_plus_in, subfolder='correlations', description=description)

    #!--------------------------------------------------------------------------------------!

    if numscalar == 1:

        print(">>> Plotting spanwise scalar auto-correlation coefficient Css(rz).")
        print(">>> Folder: plots/correlations/.")
        print()

        # Cssz
        fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

        # Limits for axes
        xliminf = 0.0
        xlimsup = Lz_plus / 2.0
        yliminf = np.min(Rssz[c,:])*1.2
        ylimsup = np.max(Rssz[c,:])*1.2

        # Auto-correlation coefficient for phi' (s: scalar)
        ax.scatter(rz, Rssz[c,:nz], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

        # Plot horizontal line at Css = 0
        ax.hlines(y=0.0, xmin=xliminf, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed')
                   
        # Axes labels
        ax.set_xlabel(r'$r_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
        ax.set_ylabel(r'$C_{\varphi \varphi}(r_z^+)$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

        # Set the plot parameters using the function 'set_plot_settings'
        # Last argument is the switcher for semilog plot (1: yes, 0: no)
        set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)
        
        # Description of .pdf file
        description = 'Auto-correlation coefficient for scalar field in spanwise direction.'

        # Save and show the figure
        save_and_show_plot('Cssz', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, y_plus_in=y_plus_in, subfolder='correlations', description=description)

    #!--------------------------------------------------------------------------------------!

#!--- Turbulent Kinetic Energy (TKE) budgets ---!
if post_tke_eq:

    """ 
    Ratio between turbulent production and dissipation;
    resizing of arrays.
    
    """
    p_eps_ratio_tke = np.zeros(ny)

    tke_turbt  = tke_turbt [:ny]
    tke_presst = tke_presst[:ny]
    tke_difft  = tke_difft [:ny]
    tke_prod   = tke_prod  [:ny]
    tke_pseps  = tke_pseps [:ny]
    
    p_eps_ratio_tke = np.divide(tke_prod,tke_pseps)

    #!--------------------------------------------------------------------------------------!

    print(">>> Plotting ratio of production over dissipation of TKE.")
    print(">>> Folder: plots/tke_stats/.")
    print()

    # Ratio of production over dissipation of TKE
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

    # Limits for axes
    xliminf = 0.1
    yliminf = 0.0
    ylimsup = 2.0

    # P/eps of TKE
    ax.scatter(y_plus[:ny], p_eps_ratio_tke[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
    # Description of .pdf file
    description = 'Ratio of production over dissipation of Turbulent Kinetic Energy (TKE). Reference data '
    
    # TTBL
    if itype == 13:

        xlimsup = 520.0
        
        # Lee & Moser (2015)
        ax.plot(y_plus_lm1000, p_eps_ratio_lm1000 + 1.0, color='C1', linestyle='-', linewidth=pp.lw)
        
        # Completing description
        description += 'Lee & Moser (2015).'
                
    # Channel    
    elif itype == 3:

        xlimsup = 300.0
        
        # Moser et al. (1999)
        ax.plot(y_plus_moser_1999, p_eps_ratio_moser_1999, color='C1', linestyle='-', linewidth=pp.lw)
        
        # Completing description
        description += 'Moser et al. (1999).'
            
    # Axes labels
    ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$P/\varepsilon$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

    # Save and show the figure
    save_and_show_plot('p_eps_ratio_tke', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, subfolder='tke_stats', description=description)
    
    #!--------------------------------------------------------------------------------------!
    
    print(">>> Plotting TKE budget terms.")
    print(">>> Folder: plots/tke_stats/.")
    print()

    # TKE budget terms
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

    # Limits for axes
    xliminf =  0.1
    xlimsup = y_plus[ny-1]*1.5
    yliminf = -0.3
    ylimsup =  0.4
    
    # Description of .pdf file
    description = 'Budget terms for Turbulent Kinetic Energy (TKE) equation. Reference data Mansour et al. (1988).'

    # Temporary: ask the user a rescale factor for tke_difft
    scale_factor = np.float64(input(">>> Scale factor for TKE diffusion: "))

    tke_difft = tke_difft*scale_factor

    # Transport terms
    ax.scatter(y_plus[:ny], -tke_turbt [:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    ax.scatter(y_plus[:ny], -tke_presst[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C1')
    ax.scatter(y_plus[:ny],  tke_difft [:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C2')
    ax.scatter(y_plus[:ny],  tke_prod  [:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C3')
    ax.scatter(y_plus[:ny], -tke_pseps [:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C4')

    # Mansour et al. (1988)            
    ax.plot(y_plus_tke_turbt_mansour,   tke_turbt_mansour,  color='C0', linestyle='-', linewidth=pp.lw)
    ax.plot(y_plus_tke_presst_mansour,  tke_presst_mansour, color='C1', linestyle='-', linewidth=pp.lw)
    ax.plot(y_plus_tke_difft_mansour,   tke_difft_mansour,  color='C2', linestyle='-', linewidth=pp.lw)
    ax.plot(y_plus_tke_prod_mansour,    tke_prod_mansour,   color='C3', linestyle='-', linewidth=pp.lw)
    ax.plot(y_plus_tke_pseps_mansour,   tke_pseps_mansour,  color='C4', linestyle='-', linewidth=pp.lw)
                   
    # Axes labels
    ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

    # Save and show the figure
    save_and_show_plot('tke_budget', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau, subfolder='tke_stats', description=description)
    
    #!--------------------------------------------------------------------------------------!



