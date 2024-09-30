
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this script, we perform:                              
!              - plotting of streamwise friction coefficient vs time       
!                and its mean value calculation (Channel);                
!              - plotting of streamwise friction coefficient vs time,     
!                friction Reynolds number vs time and streamwise          
!                friction coefficient vs friction Reynolds number (TTBL);
!              - calculation with an external subroutine of 6th order 
!                integrals of TTBL thickness parameters (delta*, theta).
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '../../4-common', 'python_common'))
sys.path.append(config_path)

# Import the plotting_params module
import plot_params as pp

# Import functions to setting up, save and show plots a
from plot_subs import set_plot_settings, save_and_show_plot

# Import function to read 'input.i3d' and 'post.prm' files and reference data
from read_files import read_input_files, read_ref_data_temp_evol

# Import function to setup flow parameters 
from set_flow_parameters import set_flow_parameters

# Import function to calculate TTBL thickness parameters at the 6th order
from cf_monitoring_subs import calculate_thickness_param

#!--------------------------------------------------------------------------------------!

# Print to screen what the program does
print()
print("!--- 'cf_monitoring.py' ---!")
print()
print(" Calculation and plotting of:")
print()
print(" Channel: ")
print(" - streamwise friction coefficient vs time ")
print("   and its mean value calculation.         ")
print()
print(" TTBL: ")
print(" - streamwise friction coefficient vs time;")
print(" - friction Reynolds number vs time;       ")
print(" - streamwise friction coefficient vs friction Reynolds number;")
print(" - 6th order accurate integrals of TTBL thickness parameters (delta*, theta).")                
print()

#!--------------------------------------------------------------------------------------!

# Create folders to store later results (e.g. cf_mean and plot)
os.makedirs('data_post',               mode=0o777, exist_ok=True)
os.makedirs('data_post/cf_monitoring', mode=0o777, exist_ok=True)
os.makedirs('plots',                   mode=0o777, exist_ok=True)
os.makedirs('plots/time_evolution',    mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
(itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, ioutput, ioutput_cf, iswitch_wo,  
 add_string, file1, filen, icrfile, nr, post_mean, post_vort, post_diss, post_corz, post_tke_eq
) = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!

#!--- Read reference data ---!
(t_gboga, 
 retau_vs_time_gboga,
 retau_gboga,
 retheta_gboga,
 utau_gboga,
 delta_99_gboga,
 disp_t_gboga,
 cf_gboga    
 ) = read_ref_data_temp_evol()

#!--- Parameters ---!
uwall, nu, twd = set_flow_parameters(itype, re)
                         
#!--- Reading of files section ---!

# Reading of grid points
y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)

# Channel (only 1 flow realization generally)
if itype == 3:
   
    # Read cf data from /data folder
    M = np.loadtxt('data/monitoring/cf_history.txt', skiprows=1, delimiter=',', dtype=np.float64)

    # Extracting quantities from the full matrix
    cfx       = M[:,3] 
    time_unit = M[:,6]

# TTBL
elif itype == 13:

    print()
    print(">>> Averaging 'cf_history.txt' files over different flow realizations.")
    
    # Do loop over different realizations
    for i in range(1, nr + 1, 1):

        # Read cf data from /data_ri folder
        M = np.loadtxt(f'data_r{i:01d}/monitoring/cf_history.txt', skiprows=1, delimiter=',', dtype=np.float64)
  
        # Extracting quantities from the full matrix
        sh_veltot = M[:,0]     # total shear velocity
        sh_velx   = M[:,1]     # streamwise shear velocity
        a_fact    = M[:,4]     # Reynolds analogy factor
        time_unit = M[:,6]     # time unit
        ts        = M[:,7]     # time step
        power_in  = M[:,10]    # power input
                             
        """
        Initialize arrays for sum
        Shear velocities, sq_sum: sum of the square, gradient multiplied by kinematic viscosity.
        This is made in order to recover the exact arithmetic average 
        (that is between gradients) (summation between linear quantities).
        """
        if i == 1:

            sh_veltotsq_sum = np.zeros(len(time_unit))    
            sh_velxsq_sum   = np.zeros(len(time_unit))
            a_fact_sum      = np.zeros(len(time_unit))
            power_in_sum    = np.zeros(len(time_unit))

        # Average the square of the total shear velocity over the realizations 
        sh_veltotsq_sum = sh_veltotsq_sum + (sh_veltot**2 / nr)

        # Average the square of the longitudinal shear velocity over the realizations 
        sh_velxsq_sum = sh_velxsq_sum + (sh_velx**2 / nr)

        # Average the Reynolds analogy factor over the realizations
        a_fact_sum = a_fact_sum + a_fact / nr
                
        # Average the power input over the realizations
        # (to be refined with the explicit calculation with wall velocities and shear velocities)
        power_in_sum = power_in_sum + power_in / nr 
                
    # Finalize the averages
    sh_veltot = np.sqrt(sh_veltotsq_sum)
    sh_velx   = np.sqrt(sh_velxsq_sum)
    a_fact    = a_fact_sum
    power_in  = power_in_sum
        
    # Calculate (streamwise) friction coefficient
    cfx = 2.0 * (sh_velx / uwall)**2
    
    #!--------------------------------------------------------------------------------------------------------!
    
    print(">>> Saving 'cf_history_realiz.txt' in data_post/cf_monitoring/.")
    print()

    # Create the file and write  
    with open('data_post/cf_monitoring/cf_history_realiz.txt', 'w') as f:
        f.write(f"{'sh_veltot (O(6))':>{pp.c_w}}, " +
                f"{'sh_velx (O(6))':>{pp.c_w}}, "   +
                f"{'cfx (O(6))':>{pp.c_w}}, "       +
                f"{'P_in':>{pp.c_w}}, "             +
                f"{'A_fact':>{pp.c_w}}, "           +
                f"{'time_unit':>{pp.c_w}}\n"        )

        for j in range(0, len(time_unit)):
            f.write(f"{sh_veltot[j]:{pp.fs6}}, "    +
                    f"{sh_velx[j]:{pp.fs6}}, "      +    
                    f"{cfx[j]:{pp.fs8}}, "          +
                    f"{power_in[j]:{pp.fs6}}, "     +
                    f"{a_fact[j]:{pp.fs}}, "        +
                    f"{time_unit[j]:{pp.fs}}\n"     )

    #!--- Section on check of mesh spacings ---!

    print()
    print(">>> Calculating grid spacings and viscous time scale at maximum cf.")

    """
    Maximum (total) shear velocity.
    Excluding first 5 savings to avoid the IC peak of cf.
    The same index is used to plot, in order to match what we plot
    to what we use to find the cf maximum.
    """
    
    # Lower index cf (to avoid first points in plotting and when we find the maximum) 
    # (l: lower; i: index)
    # First 5 points avoided in this manner.
    li_cf = 4

    max_sh_veltot = np.max(sh_veltot[li_cf:]) 
    
    # Related viscous lengths
    delta_nu_tot = nu / max_sh_veltot

    # Mesh spacings (dimensional)
    delta_x  = Lx / nx
    delta_yw = y[1]  
    delta_z  = Lz / nz

    """
    Non-dimensional mesh spacings and viscous time unit
    (p: plus, tot: total shear velocity).
    """
    delta_x_p_tot  = delta_x  / delta_nu_tot
    delta_yw_p_tot = delta_yw / delta_nu_tot
    delta_z_p_tot  = delta_z  / delta_nu_tot

    t_nu_tot = nu / max_sh_veltot**2
    
    # Minimum viscous time unit
    t_nu_min = np.min(t_nu_tot)
    
    # Ratio of time-step dt and minimum viscous time unit
    dt_plus = round(dt / t_nu_min, 4)

    # Saving "num_resolutions.txt" (num: numerical)
    print(">>> Saving 'num_resolutions.txt' in data_post/cf_monitoring/.")
    print()

    # Write and save to .txt file 
    with open('data_post/cf_monitoring/num_resolutions.txt', 'w') as f:
        f.write('Maximum non-dimensional grid spacings and minimum viscous time scale.\n')
        f.write('Rescaling with total shear velocity,\n')
        f.write('based on the norm of the wall shear stress vector.\n')
        f.write('\n')
        f.write(f'Time-step dt in viscous units, dt^+ = dt / t_nu_min = {dt_plus}.\n')
        f.write('\n')
        f.write(f"{'delta_x+_max':>{pp.c_w}}, "  +
                f"{'delta_yw+_max':>{pp.c_w}}, " +
                f"{'delta_z+_max':>{pp.c_w}}, "  + 
                f"{'t_nu_max':>{pp.c_w}}\n"      )
        
        f.write(f"{delta_x_p_tot:{pp.fs6}}, "    +
                f"{delta_yw_p_tot:{pp.fs6}}, "   +
                f"{delta_z_p_tot:{pp.fs6}}, "    +
                f"{t_nu_tot:{pp.fs6}}\n"         )

    """
    Call subroutine for calculations of 6th order TTBL thickness parameters,
    mean_stats runtime averaged with different flow realizations and 
    TTBL thickness delta_99.     
    """
    (delta_99, disp_t, mom_t) = calculate_thickness_param(sh_veltot,sh_velx)
    
    # Calculate the (streamwise) friction Reynolds number (averaged over the realizations)
    re_tau = sh_velx * delta_99 / nu

    # Calculate Reynolds numbers based on displacement thickness and momentum thickness
    re_disp_t = disp_t * re
    re_mom_t  = mom_t  * re

    print(">>> Saving 'high_order_integrals_evol.txt' in data_post/cf_monitoring/.")
    print()

    # Create the file and write  
    with open('data_post/cf_monitoring/high_order_integrals_evol.txt', 'w') as f:
        f.write(f"{'cfx':>{pp.c_w}}, "          +
                f"{'delta_99':>{pp.c_w}}, "     +
                f"{'disp_t':>{pp.c_w}}, "       +
                f"{'mom_t':>{pp.c_w}}, "        +  
                f"{'Re_tau':>{pp.c_w}}, "       +
                f"{'Re_delta*':>{pp.c_w}}, "    +
                f"{'Re_theta':>{pp.c_w}}, "     +
                f"{'time_unit':>{pp.c_w}}\n"    )

        for j in range(0, len(time_unit)):
            f.write(f"{cfx[j]:{pp.fs8}}, "      +
                    f"{delta_99[j]:{pp.fs}}, "  +
                    f"{disp_t[j]:{pp.fs}}, "    +
                    f"{mom_t[j]:{pp.fs}}, "     + 
                    f"{re_tau[j]:{pp.fs}}, "    +
                    f"{re_disp_t[j]:{pp.fs}}, " +
                    f"{re_mom_t[j]:{pp.fs}}, "  +
                    f"{time_unit[j]:{pp.fs}}\n" )

#!--------------------------------------------------------------------------------------!

# Only for a channel
if itype == 3:

    # Calculating the delta (in time units) between savings (skipping first saving that is IC usually)
    a = time_unit[1]
    b = time_unit[2]
    delta = b - a
    
    # Asking the user for lower limit for the range of time units for cf
    lower_tu = np.float64(input("Specify a lower range for time units (T): "))
    print()

    # Calculating its related index and show it
    lower_index = int(lower_tu / delta)
    print("Correspondent snapshot index:", lower_index)
    print()

    # Average (lower TU is included)
    mean_cf = np.mean(cfx[lower_index:])
    print("Mean cf value:", mean_cf)
    print()
    
    # Number of snapshots used and total average time (lower TU is included)
    last_index = len(time_unit) 
    n_snap = last_index - lower_index
    t_tot = (n_snap - 1)*delta

#!--------------------------------------------------------------------------------------!

print()
print(">>> Plotting streamwise friction coefficient as function of time.")
print()

# Axes ranges
xliminf = time_unit[0]
xlimsup = time_unit[-1]
yliminf = np.min(cfx) * 0.0
ylimsup = np.max(cfx) * 1.2 

# Subplots environment
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
   
# TTBL only
if itype == 13:

    # Friction coefficient
    ax.plot(time_unit, cfx, color='C0', linestyle='-', linewidth=pp.lw)

    # x-axis label
    ax.set_xlabel(r'$t U_w/D$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Channel flow only
elif itype == 3:

    # Friction coefficient
    ax.scatter(time_unit, cfx, marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

    # x-axis label 
    ax.set_xlabel(r'$t\frac{U_p}{h}$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    # Create a rectangle patch to show points we are excluding from average
    rect = patches.Rectangle((0, 0), lower_tu, ylimsup, linewidth=0, edgecolor='none', facecolor='r', alpha=0.1)
    
    # Add the patch to the plot
    ax.add_patch(rect)
    
    # Horizontal line to show mean cf value
    ax.hlines(y=mean_cf, xmin=lower_tu, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed', label=f'Mean value: {mean_cf:.3e}')
               
    # Create the file and write  
    with open('data_post/cf_monitoring/cf_mean.txt', 'w') as f:
        f.write(f"{'cf_mean':<{pp.c_w}}, "  +
                f"{'t_tot':<{pp.c_w}}, "    +
                f"{'delta_TU':<{pp.c_w}}, " +
                f"{'n_snap':<{pp.c_w}}\n"   )

        f.write(f"{mean_cf:{pp.fs8}}, "     +
                f"{t_tot:{pp.fs}}, "        +
                f"{delta:{pp.fs}}, "        +
                f"{n_snap:{pp.fs}}\n"       )

# y-axis label
ax.set_ylabel(r'$c_{f x}$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Description
description = 'Time evolution of streamwise friction coefficient.'
    
# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)
    
# Save and show the figure
save_and_show_plot('cfx_vs_time', add_string=add_string, subfolder='time_evolution', description=description)

#!--------------------------------------------------------------------------------------!

# TTBL only
if itype == 13:

    print()
    print(">>> Plotting (streamwise) friction Reynolds number as function of time.")
    print(">>> Reference data Cimarelli et al. (2024a), data with total shear velocity.")
    print()

    # Subplots environment
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
   
    # Friction Reynolds number
    ax.plot(time_unit, re_tau, color='C0', linestyle='-', linewidth=pp.lw)
    
    # G. Boga (Cimarelli et al. (2024a))
    ax.plot(t_gboga, retau_vs_time_gboga, color='C1', linestyle='-', linewidth=pp.lw)
    
    # Axes labels
    ax.set_xlabel(r'$t U_w/D$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$Re_\tau$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

    # Axes ranges
    xliminf = time_unit[0]
    xlimsup = time_unit[-1]
    yliminf = np.min(re_tau) * 0.0
    ylimsup = np.max(re_tau) * 1.2 
              
    # Description
    description  = 'Time evolution of (streamwise) friction Reynolds number. '
    description += 'Reference data Cimarelli et al. (2024a).'
    
    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)
    
    # Save and show the figure
    save_and_show_plot('retaux_vs_time', add_string=add_string, subfolder='time_evolution', description=description)
    
    #!--------------------------------------------------------------------------------------!

    print()
    print(">>> Plotting streamwise friction coefficient as function of (streamwise) friction Reynolds number.")
    print(">>> Reference data Cimarelli et al. (2024a), data with total shear velocity.")
    print()

    # Subplots environment
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
   
    # Friction coefficient (li_cf: lower index cf, to exclude the peak of IC)
    ax.plot(re_tau[li_cf:], cfx[li_cf:], color='C0', linestyle='-', linewidth=pp.lw)
    
    # G. Boga (Cimarelli et al. (2024a)) 
    ax.plot(retau_gboga, cf_gboga, color='C1', linestyle='-', linewidth=pp.lw)
    
    # Axes labels
    ax.set_xlabel(r'$Re_\tau$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$c_{f}$',   fontsize=pp.fla, labelpad=pp.pad_axes_lab)

    # Axes ranges
    xliminf = np.min(re_tau)
    xlimsup = np.max(re_tau)
    yliminf = np.min(cfx) * 0.0
    ylimsup = np.max(cfx) * 1.2 
              
    # Description
    description  = 'Time evolution of streamwise friction coefficient against '\
                   '(streamwise) friction Reynolds number. '
    description += 'Reference data Cimarelli et al. (2024a).'
    
    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)
    
    # Save and show the figure
    save_and_show_plot('cfx_vs_retaux', add_string=add_string, subfolder='time_evolution', description=description)

    #!--------------------------------------------------------------------------------------!

    print()
    print(">>> Plotting friction Reynolds number vs momentum thickness Reynolds number.")
    print(">>> Reference data Cimarelli et al. (2024a).")
    print(">>> Reference straight line Schlatter & Orlu (2010).")
    print()

    # Subplots environment
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
   
    # Friction Reynolds number vs momentum thickness Reynolds number
    ax.plot(re_mom_t, re_disp_t, color='C0', linestyle='-', linewidth=pp.lw)
    
    # G. Boga (Cimarelli et al. (2024a)) 
    ax.plot(retheta_gboga, retau_gboga, color='C1', linestyle='-', linewidth=pp.lw)

    # Schlatter & Orlu (2010), straight line (so: Schlatter & Orlu)
    re_theta_so = np.linspace(1, 3000, 3000)
    re_tau_so   = 1.13 * re_theta_so**0.843
    ax.plot(re_theta_so, re_tau_so, color='C2', linestyle='-', linewidth=pp.lw)
    
    # Axes labels
    ax.set_xlabel(r'$Re_\theta$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$Re_\tau$',   fontsize=pp.fla, labelpad=pp.pad_axes_lab)

    # Axes ranges
    xliminf = np.min(re_mom_t)
    xlimsup = np.max(re_mom_t)
    yliminf = np.min(re_tau) * 0.0
    ylimsup = np.max(re_tau) * 1.2 
              
    # Description
    description  = 'Friction Reynolds number against '\
                   'momentum thickness Reynolds number. '
    description += 'Reference data Cimarelli et al. (2024a).'
    description += 'Reference straight line Schlatter & Orlu (2010).'
    
    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)
    
    # Save and show the figure
    save_and_show_plot('re_tau_vs_re_theta', add_string=add_string, subfolder='time_evolution', description=description)

    #!--------------------------------------------------------------------------------------!

    print()
    print(">>> Plotting displacement thickness Reynolds number vs momentum thickness Reynolds number.")
    print(">>> Reference data Cimarelli et al. (2024a).")
    print()

    # Subplots environment
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
   
    # Displacement thickness Reynolds number vs momentum thickness Reynolds number
    ax.plot(re_mom_t, re_disp_t, color='C0', linestyle='-', linewidth=pp.lw)
    
    # G. Boga (Cimarelli et al. (2024a)) 
    # (displacement thickness multiplied by trip Reynolds number, that is the same in both cases).
    ax.plot(retheta_gboga, disp_t_gboga*re, color='C1', linestyle='-', linewidth=pp.lw)
    
    # Axes labels
    ax.set_xlabel(r'$Re_\theta$',     fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$Re_{\delta^*}$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

    # Axes ranges
    xliminf = np.min(re_mom_t)
    xlimsup = np.max(re_mom_t)
    yliminf = np.min(re_disp_t) * 0.0
    ylimsup = np.max(re_disp_t) * 1.2 
              
    # Description
    description  = 'Displacement thickness Reynolds number against '\
                   'momentum thickness Reynolds number. '
    description += 'Reference data Cimarelli et al. (2024a).'
    
    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)
    
    # Save and show the figure
    save_and_show_plot('re_disp_t_vs_re_theta', add_string=add_string, subfolder='time_evolution', description=description)

#!--------------------------------------------------------------------------------------!

print()
print(">>> End.")
print()



