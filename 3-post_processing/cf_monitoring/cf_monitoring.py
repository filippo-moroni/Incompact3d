#!---------------------------------------------------------!
#! With this script, we perform plotting of friction       !
#! coefficient and its mean value calculation.             !
#!---------------------------------------------------------!

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '..', 'python_common'))
sys.path.append(config_path)

# Import the plotting_params module
import plot_params as pp

# Import function to set plots
from plot_settings import set_plot_settings

# Import function to read 'input.i3d' and 'post.prm' files
from read_incompact3d_files import read_input_files

#!--------------------------------------------------------------------------------------!

# Default value for axes limits
xliminf     = 0.0              # x axis inferior limit
xlimsup     = 1500.0           # x axis superior limit
yliminf     = 0.0              # y axis inferior limit
ylimsup     = 0.008            # y axis superior limit

#!--------------------------------------------------------------------------------------!

# Create folders to store later results (e.g. cf_mean and plot)
os.makedirs('data_post', mode=0o777, exist_ok=True)
os.makedirs('plots',     mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
itype, nx, ny, nz, Lx, Ly, Lz, re, iswitch_wo, file1, filen, icrfile, nr, add_string = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!
                         
#!--- Reading of files section ---!

# Path for generic data
data_path = "/data"

# Check if the path exists and is a directory
if os.path.exists(data_path) and os.path.isdir(data_path):
    print("The folder exists.")
    
    # Read cf data from /data folder
    M1 = np.loadtxt('data/monitoring/cf_history.txt', skiprows=1, delimiter=',', dtype=np.float64)

else:
    print("The folder does not exist.")
    
    # Read cf data from /data_ri folder, where the ri (realization number) is read from 'post.prm'
    M1 = np.loadtxt(f'data_r{nr:01d}/monitoring/cf_history.txt', skiprows=1, delimiter=',', dtype=np.float64)
  
# Extracting quantities from the full matrix
cfx       = M1[:,4] 
time_unit = M1[:,7] 

# Axes ranges
xliminf = time_unit[0]
xlimsup = time_unit[-1]
yliminf = np.min(cfx) * 0.0
ylimsup = np.max(cfx) * 1.2

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

#!--- Plot section, friction coefficient ---!
print()
print("!--- Plotting of friction coefficient ---!")
print()

# Subplots environment
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
   
# Friction coefficient
ax.scatter(time_unit, cfx, marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

# Axes labels
ax.set_xlabel(r'$t$',       fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$c_{f,x}$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Channel flow only
if itype == 3:

    # x-axis label 
    ax.set_xlabel(r'$t\frac{U_p}{h}$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    # Create a rectangle patch to show points we are excluding from average
    rect = patches.Rectangle((0, 0), lower_tu, ylimsup, linewidth=0, edgecolor='none', facecolor='r', alpha=0.1)
    
    # Add the patch to the plot
    ax.add_patch(rect)
    
    # Horizontal line to show mean cf value
    ax.hlines(y=mean_cf, xmin=lower_tu, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed', label=f'Mean value: {mean_cf:.3e}')
               
    # Create the file and write  
    with open('data_post/cf_mean.txt', 'w') as f:
        f.write(f"{'cf_mean':<{pp.c_w}}, "  +
                f"{'t_tot':<{pp.c_w}}, "    +
                f"{'delta_TU':<{pp.c_w}}, " +
                f"{'n_snap':<{pp.c_w}}\n"   )

        f.write(f"{mean_cf:{pp.fs2}}, "     +
                f"{t_tot:{pp.fs}}, "        +
                f"{delta:{pp.fs}}, "        +
                f"{n_snap:{pp.fs}}\n"       )
               
               
# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Saving the figure and show it
plt.savefig(f'plots/cf_vs_time_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)
plt.show()

#!--------------------------------------------------------------------------------------!



