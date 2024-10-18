
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this script, we perform:
!              - averages of runtime mean statistics for a TTBL                               
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

# Import function to average runtime mean stats with different flow realizations
from cf_monitoring_subs import average_runtime_mean_stats

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
print(" - 6th order accurate integrals of TTBL thickness parameters (delta*, theta);")
print(" - friction Reynolds number as function of momentum thickness Reynolds number;")
print(" - displacement thickness Reynolds number as function of momentum thickness Reynolds number.")
print()

#!--------------------------------------------------------------------------------------!

# Create folders to store later results
os.makedirs('time_evolution',       mode=0o777, exist_ok=True)
os.makedirs('num_resolutions',      mode=0o777, exist_ok=True)
os.makedirs('plots',                mode=0o777, exist_ok=True)
os.makedirs('plots/time_evolution', mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
(itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, itimescheme, ioutput, ioutput_cf, iswitch_wo,  
 add_string, file1, filen, icrfile, nr, post_mean, post_grad, post_corz, post_tke_eq
) = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!

#!--- Read reference data ---!
(t_gboga, retau_vs_time_gboga, retau_gboga, retheta_gboga, utau_gboga, delta_99_gboga, disp_t_gboga, cf_gboga    
) = read_ref_data_temp_evol()

#!--- Parameters and mesh ---!
(uwall, nu, twd, y) = set_flow_parameters(re)
                         
# Start of calculations for a TTBL
if itype == 13:

    print()
    print(">>> Opening 'cf_history.txt' files over different flow realizations.")
    print()
        
    # Number of savings due to 'print_cf' subroutine of Incompact3d solver 'modified';
    # We can restrict the range to plot in order to be more flexible (e.g. if simulations have not been completed).
    # 'ilast' can be changed from the input file 'input.i3d'.
    nsavings = ilast // ioutput_cf + 1
    
    # Asking to the user the width of the time-window average (number of indexes)
    print(">>> Enter the number of snapshots to be used for centered time-window average: ")
    print("    0   : time-window average is not performed;")
    print("    twi : time-window index (half number of snapshots excluding the central one).")
    print()
    
    # Input from the user
    time_window_index = int(input("Time window index = "))
    
    # Take alias
    twi = time_window_index
            
    # Number of snapshots used to window averaging in time
    nt = (2 * time_window_index) + 1
    
    # Denominator of the divisions (number of snapshots in the time window average times number of flow realizations) 
    den = nt * nr
    
    # Number of elements in the reduced arrays (number of savings reduced) (red: reduced)
    nsavings_red = nsavings - 2*time_window_index
        
    # Initialise arrays for sum
    sh_vel_tot_sq_sum = np.zeros(nsavings_red)    
    sh_vel_x_sq_sum   = np.zeros(nsavings_red)
    mg_phi_w_sum      = np.zeros(nsavings_red)
    a_fact_sum        = np.zeros(nsavings_red)
    time_unit_sum     = np.zeros(nsavings_red)
    power_in_sum      = np.zeros(nsavings_red)
      
    # Do loop over different realizations
    for i in range(1, nr + 1, 1):

        # Read cf data from /data_ri folder
        M = np.loadtxt(f'data_r{i:01d}/monitoring/cf_history.txt', skiprows=1, delimiter=',', dtype=np.float64)
  
        # Extracting quantities from the full matrix
        sh_vel_tot = M[:nsavings,0]     # total shear velocity
        sh_vel_x   = M[:nsavings,1]     # streamwise shear velocity
        mg_phi_w   = M[:nsavings,4]     # mean scalar gradient at the wall
        a_fact     = M[:nsavings,5]     # Reynolds analogy factor
        time_unit  = M[:nsavings,7]     # time unit
        ts         = M[:nsavings,8]     # time step
        power_in   = M[:nsavings,11]    # power input
                             
        """
        !--------------------------------------------------------------------------------------------------!
         For shear velocities, we sum their square value, in order to recover the exact arithmetic average. 
         (Summation between linear quantities, gradients multiplied by kinematic viscosity).
        !--------------------------------------------------------------------------------------------------!
        
        """
        
        # Cycle on ti: time index, with restriction at the limits of the interval if time-window averaging is enabled
        for ti in range(0+twi, nsavings-twi, 1):
            
            # Time-window average cycle (we are summing centered on time index 'ti')
            for i in range(-twi, twi+1, 1):
                                        
                # Average the square of the total shear velocity over the realizations 
                sh_vel_tot_sq_sum[ti-twi] = sh_vel_tot_sq_sum[ti-twi] + sh_vel_tot[ti+i]**2 / den
                
                # Average the square of the longitudinal shear velocity over the realizations 
                sh_vel_x_sq_sum[ti-twi] = sh_vel_x_sq_sum[ti-twi] + sh_vel_x[ti+i]**2 / den
                
                # Average the mean scalar gradient at the wall
                mg_phi_w_sum[ti-twi] = mg_phi_w_sum[ti-twi] + mg_phi_w[ti+i] / den
                
                # Average the Reynolds analogy factor over the realizations
                a_fact_sum[ti-twi] = a_fact_sum[ti-twi] + a_fact[ti+i] / den
                
                # Time-unit sum in order to check time-window averaging
                time_unit_sum[ti-twi] = time_unit_sum[ti-twi] + time_unit[ti+i] / den
                
                # Average the power input over the realizations
                power_in_sum[ti-twi] = power_in_sum[ti-twi] + power_in[ti+i] / den 
    
    # Resize arrays
    sh_vel_tot_sq = np.zeros(nsavings_red)    
    sh_vel_x_sq   = np.zeros(nsavings_red)
    mg_phi_w      = np.zeros(nsavings_red)
    a_fact        = np.zeros(nsavings_red)
    time_unit     = np.zeros(nsavings_red)
    power_in      = np.zeros(nsavings_red)
                
    # Finalize the averages
    sh_vel_tot = np.sqrt(sh_vel_tot_sq_sum)
    sh_vel_x   = np.sqrt(sh_vel_x_sq_sum)
    mg_phi_w   = mg_phi_w_sum
    time_unit  = time_unit_sum
    a_fact     = a_fact_sum
    power_in   = power_in_sum
        
    # Calculate (streamwise) friction coefficient
    cfx = 2.0 * (sh_vel_x / uwall)**2
    
    print()
    print(">>> Average of shear velocities, mean scalar gradient, Reynolds analogy factor and power input: done.")
    
    #!--------------------------------------------------------------------------------------------------------!
    
    """
    !---------------------------------------------------------------------------------------------------------!
    Call subroutine for the calculations of:
     - mean_stats runtime averaged with different flow realizations (and time-window average possibility);
     - 6th order integrals for TTBL thickness parameters delta* and theta;
     - TTBL thickness delta_99;
     - maximum mesh spacing in y-direction at the BL interface in viscous units.
    !---------------------------------------------------------------------------------------------------------!     
    """
    (delta_99, disp_t, mom_t, max_delta_yd_plus) = average_runtime_mean_stats(sh_vel_tot, sh_vel_x, mg_phi_w, nsavings, time_window_index, den)
    
    print()
    print(">>> Average of runtime mean statistics with different flow realizations")
    print("    and calculation of time evolution of TTBL thickness parameters: done.")
        
    # Calculate the (streamwise) friction Reynolds number (averaged over the realizations)
    re_tau = np.zeros(nsavings_red) 
    re_tau = sh_vel_x * delta_99 / nu

    # Calculate Reynolds numbers based on displacement thickness and momentum thickness
    re_disp_t = disp_t * re
    re_mom_t  = mom_t  * re
    
    # Calculate viscous time unit
    t_nu = np.zeros(nsavings_red)
    t_nu = nu / (sh_vel_tot**2)

    print()
    print(">>> Saving 'time_evolution.txt' in time_evolution/.")
    print("    This file stores the following quantities:")
    print("     - cfx;")
    print("     - delta_99, delta*, theta;")
    print("     - Re_tau, Re_delta*, Re_theta;")
    print("     - power input, P_in;")
    print("     - Reynolds analogy factor, A_fact.")
    print()
    print(">>> Statistics averaged on a time window of amplitude")
    print(f"   T_avg = {2.0*twi*dt*ioutput_cf},")        
    print(f"   with number of snapshots in time, nt = {nt}.") 
    print()

    # Create the file and write  
    with open('time_evolution/time_evolution.txt', 'w') as f:
        f.write('This file stores the main quantities of interest to analyse the behaviour of TTBLs.\n')
        f.write('Acronyms & quantities:\n')
        f.write(f' - cfx       : streamwise friction coefficient;\n')
        f.write(f' - delta_99  : TTBL thickness delta_99;\n')
        f.write(f' - disp_t    : TTBL displacement thickness (delta*);\n')        
        f.write(f' - mom_t     : TTBL momentum thickness (theta);\n')        
        f.write(f' - Re_tau    : friction Reynolds number;\n')        
        f.write(f' - Re_delta* : Reynolds number based on displacement thickness;\n')        
        f.write(f' - Re_theta  : Reynolds number based on momentum thickness;\n')        
        f.write(f' - P_in      : power input (valid for both fixed and oscillating wall);\n')        
        f.write(f' - A_fact    : Reynolds analogy factor.\n')
        f.write(f' - t_nu      : viscous time unit, based on total shear velocity.\n')        
        f.write(f' - time_unit : non-dimensional outer time scale, based on wall velocity and trip wire diameter.\n')            
        f.write('\n')
        f.write(f'Statistics averaged on a time window of amplitude, T_avg = {2.0*twi*dt*ioutput_cf},\n')        
        f.write(f'with number of snapshots in time, nt = {nt}.\n')                
        f.write('\n')
        f.write(f"{'cfx':>{pp.c_w}}, "          +
                f"{'delta_99':>{pp.c_w}}, "     +
                f"{'disp_t':>{pp.c_w}}, "       +
                f"{'mom_t':>{pp.c_w}}, "        +  
                f"{'Re_tau':>{pp.c_w}}, "       +
                f"{'Re_delta*':>{pp.c_w}}, "    +
                f"{'Re_theta':>{pp.c_w}}, "     +
                f"{'P_in':>{pp.c_w}}, "         +
                f"{'A_fact':>{pp.c_w}}, "       +
                f"{'t_nu':>{pp.c_w}}, "         +                
                f"{'time_unit':>{pp.c_w}}\n"    )

        for j in range(0, nsavings_red):
            f.write(f"{cfx[j]:{pp.fs8}}, "      +
                    f"{delta_99[j]:{pp.fs}}, "  +
                    f"{disp_t[j]:{pp.fs}}, "    +
                    f"{mom_t[j]:{pp.fs}}, "     + 
                    f"{re_tau[j]:{pp.fs}}, "    +
                    f"{re_disp_t[j]:{pp.fs}}, " +
                    f"{re_mom_t[j]:{pp.fs}}, "  +
                    f"{power_in[j]:{pp.fs6}}, " +
                    f"{a_fact[j]:{pp.fs}}, "    +
                    f"{t_nu[j]:{pp.fs}}, "      +                    
                    f"{time_unit[j]:{pp.fs}}\n" )
    
    #!--------------------------------------------------------------------------------------------------------!
    
    #!--- Section on check of mesh spacings ---!

    print()
    print(">>> Calculating maximum grid spacings in wall units and minimum viscous time scale.")
    print(">>> We expect to observe maximum grid spacings along x and z and minimum viscous time at peak cf,")
    print("    while maximum wall-normal grid spacing at the BL interface at the end of the simulation.")

    # Maximum (total) shear velocity
    max_sh_vel_tot = np.max(sh_vel_tot) 
    
    # Related viscous length
    delta_nu_tot = nu / max_sh_vel_tot

    # Mesh spacings (dimensional)
    delta_x  = Lx / nx
    delta_yw = y[1]  
    delta_z  = Lz / nz

    """
    !----------------------------------------------------!
     Non-dimensional mesh spacings and viscous time unit
     (p: plus, tot: total shear velocity).
    !----------------------------------------------------!
    
    """
            
    delta_x_p_tot  = delta_x  / delta_nu_tot
    delta_yw_p_tot = delta_yw / delta_nu_tot
    delta_z_p_tot  = delta_z  / delta_nu_tot

    t_nu_min_tot = nu / max_sh_vel_tot**2
        
    # Ratio of time-step dt and minimum viscous time unit
    dt_plus = round(dt / t_nu_min_tot, 4)

    # Saving "num_resolutions.txt" (num: numerical)
    print(">>> Saving 'worst_num_resolutions.txt' in num_resolutions/.")
    print()

    # Write and save to .txt file 
    with open('num_resolutions/worst_num_resolutions.txt', 'w') as f:
        f.write('Maximum non-dimensional grid spacings and minimum viscous time scale.\n')
        f.write('Rescaling with total shear velocity, based on the norm of the total wall shear stress vector.\n')
        f.write('\n')
        f.write(f'Flowcase: {add_string}.\n')
        f.write('\n')
        f.write(f'Time-step dt = {dt}.\n')        
        f.write(f'Time-step dt in viscous units, dt^+ = dt / t_nu_min = {dt_plus}.\n')
        f.write('\n')
        f.write('Abbreviations:\n')
        f.write(' - x     : streamwise direction;\n')
        f.write(' - y     : wall-normal direction;\n')
        f.write(' - z     : spanwise direction;\n')
        f.write(' - delta : mesh spacing;\n')
        f.write(' - w     : first element at the wall;\n')        
        f.write(' - d     : boundary layer interface (d: small letter greek delta);\n')
        f.write('\n')
        f.write(f"{'delta_x+_max':>{pp.c_w}}, "   +
                f"{'delta_yw+_max':>{pp.c_w}}, "  +
                f"{'delta_z+_max':>{pp.c_w}}, "   +
                f"{'delta_yd+_max':>{pp.c_w}}, "  +
                f"{'t_nu_min':>{pp.c_w}}\n"       )
        
        f.write(f"{delta_x_p_tot:{pp.fs6}}, "     +
                f"{delta_yw_p_tot:{pp.fs6}}, "    +
                f"{delta_z_p_tot:{pp.fs6}}, "     +
                f"{max_delta_yd_plus:{pp.fs6}}, " +                 
                f"{t_nu_min_tot:{pp.fs6}}\n"      )

#!--------------------------------------------------------------------------------------!

"""
!--------------------------------------------------------------------------------------!
 Calculations for a channel (at the moment, we are limited to 1 flow realization only;
 however, it is not so common to make different channel flow realizations.
!--------------------------------------------------------------------------------------!
"""

# Channel
if itype == 3:
   
    # Read cf data from /data folder
    M = np.loadtxt('data/monitoring/cf_history.txt', skiprows=1, delimiter=',', dtype=np.float64)

    # Extracting quantities from the full matrix
    cfx       = M[:,3] 
    time_unit = M[:,6]

    # Calculating the delta (in time units) between savings (skipping first saving that is IC usually)
    a = time_unit[1]
    b = time_unit[2]
    delta = b - a
    
    # Asking the user for lower limit for the range of time units for cf
    lower_tu = np.float64(input("Specify a lower range for time units (T): "))
    print()

    # Calculating its related index and show it
    lower_index = int(lower_tu / delta)
    print(">>> Correspondent snapshot index on 'cf_history.txt' file:", lower_index)
    print()

    # Average (lower TU is included) (a futher check on precision must be made, the plotted line seems to low wrt to data points)
    mean_cf = np.mean(cfx[lower_index:], dtype=np.float128)
    
    # Define 1000.0 with quad precision
    onethousand = np.float128(1000.0)
   
    # Rescale friction coefficient by a factor of 1000.0
    mean_cf = mean_cf*onethousand
    print(f">>> Mean cf value, 10^3 <cf> = {mean_cf:.3f}")
    print()
    
    # Number of snapshots used and total average time (lower TU is included)
    last_index = len(time_unit) 
    n_snap = last_index - lower_index
    t_tot = (n_snap - 1)*delta

#!--------------------------------------------------------------------------------------!

print()
print(">>> Plotting friction coefficient as function of time.")
print()

# Axes ranges
xliminf = np.min(time_unit)
xlimsup = np.max(time_unit)
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
            
    # Vertical line to show the lower time-unit we have selected
    ax.vlines(x=lower_tu, ymin=yliminf, ymax=ylimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed')

    # Text to show in the plot the lower time-unit
    ax.text(lower_tu*1.1, mean_cf/2000.0, fr'$t = {lower_tu}$', color='k', fontsize=4, ha='left')
    
    # Horizontal line to show mean cf value
    ax.hlines(y=mean_cf/onethousand, xmin=lower_tu, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed', label=f'Mean value: {mean_cf:.3e}')
    
    # Create folder to store cf_mean 
    os.makedirs('cf_stats', mode=0o777, exist_ok=True)
           
    # Create the file and write  
    with open('cf_stats/cf_mean.txt', 'w') as f:
        f.write('Average of friction coefficient for a Channel.\n')
        f.write('\n')
        f.write(f'Flowcase: {add_string}.\n')
        f.write('\n')
        f.write('Abbreviations:\n')
        f.write(' - 10^3 <cf> : friction coefficient average (times 10^3);\n')
        f.write(' - t_tot     : total time of average (outer time, based on channel half-height and centerline velocity of the related laminar Poiseuille flow);\n')
        f.write(' - delta_TU  : delta of time units (TU) between different savings of cf;\n')
        f.write(' - n_snap    : number of snapshots used in the average.\n')
        f.write('\n')
        f.write(f"{'10^3 <cf>':>{pp.c_w}}, " +
                f"{'t_tot':>{pp.c_w}}, "     +
                f"{'delta_TU':>{pp.c_w}}, "  +
                f"{'n_snap':>{pp.c_w}}\n"    )

        f.write(f"{mean_cf:{pp.fs}}, "      +
                f"{t_tot:{pp.fs}}, "        +
                f"{delta:{pp.fs}}, "        +
                f"{n_snap:{pp.fs}}\n"       )

# y-axis label
ax.set_ylabel(r'$c_f$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Description
description = 'Time evolution of (streamwise) friction coefficient.'
    
# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)
    
# Save and show the figure
save_and_show_plot('cf_vs_time', add_string=add_string, subfolder='time_evolution', description=description)

#!--------------------------------------------------------------------------------------!

# TTBL only
if itype == 13:

    print()
    print(">>> Plotting friction Reynolds number as function of time.")
    print("    Reference data Cimarelli et al. (2024a), data with total shear velocity.")
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
    save_and_show_plot('retau_vs_time', add_string=add_string, subfolder='time_evolution', description=description)
    
    #!--------------------------------------------------------------------------------------!

    print()
    print(">>> Plotting friction coefficient as function of friction Reynolds number.")
    print("    Reference data Cimarelli et al. (2024a), data with total shear velocity.")
    print()

    # Subplots environment
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
   
    # Friction coefficient
    ax.plot(re_tau, cfx, color='C0', linestyle='-', linewidth=pp.lw)
    
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
    save_and_show_plot('cf_vs_retau', add_string=add_string, subfolder='time_evolution', description=description)

    #!--------------------------------------------------------------------------------------!

    print()
    print(">>> Plotting friction Reynolds number as function of momentum thickness Reynolds number.")
    print("    Reference data Cimarelli et al. (2024a).")
    print("    Reference straight line Schlatter & Orlu (2010).")
    print()

    # Subplots environment
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
   
    # Friction Reynolds number vs momentum thickness Reynolds number
    ax.plot(re_mom_t, re_tau, color='C0', linestyle='-', linewidth=pp.lw)
    
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
    print(">>> Plotting displacement thickness Reynolds number as function of momentum thickness Reynolds number.")
    print("    Reference data Cimarelli et al. (2024a).")
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



