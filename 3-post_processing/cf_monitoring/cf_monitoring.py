#!---------------------------------------------------------!
#! With this script, we perform plotting of friction       !
#! coefficient and its mean value calculation.             !
#!---------------------------------------------------------!

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from matplotlib.ticker import (AutoMinorLocator)
import matplotlib.patches as patches

# Settings
np.seterr(divide='ignore', invalid='ignore')

plt.rcParams.update({ 
    "text.usetex": True,  
    "font.family": "serif",
    "font.sans-serif": "Computer Modern",
    "figure.autolayout": True,
})

# Parameters for plotting
lw           = 0.6             # linewidth for plots
markersize   = 8.0             # marker size for scatter plot
fla          = 10              # fontsize of labels of x and y axes (major labels, variables)
fla2         = 4.5             # fontsize of numbers of x and y axes 
pad_axes_lab = 2               # padding of axes labels
pad_numbers  = 3               # padding of numbers on both axes
lmajt        = 4               # length of major ticks
lmint        = 2               # length of minor ticks
tick_width   = 0.5             # width of ticks and external box

# Page settings (A4 paper format: 8.3 x 11.7 inches)
xinches      = 2.6             # size in inches in x direction of the image
yinches      = 2.2             # size in inches in y direction of the image

# Default value for axes limits
xliminf     = 0.0             # x axis inferior limit
xlimsup     = 1500.0          # x axis superior limit
yliminf     = 0.0             # y axis inferior limit
ylimsup     = 0.008           # y axis superior limit

# Axes width
mpl.rcParams['axes.linewidth'] = tick_width

# Set some useful colors
grey = [0.5, 0.5, 0.5]

# CPG option when CFR is imposed
cpg_check = 'F'

#!--------------------------------------------------------------------------------------!

# Additional string at the end of the filename 
add_string = input("Specify the suffix for the filename: ")
print()

# Read if we are plotting a channel or a TTBL
with open('input.i3d', 'r') as file:
    
    # Read all lines into a list
    lines = file.readlines()
    
    # Extract the 8th line, where itype is specified 
    line = lines[7]  
    
    # Removing characters in front of the itype value
    itype = line.split('=')[-1].strip()
    
    # Convert to integer
    itype = int(itype)

# If channel flow, see if CPG option is enabled or not
if itype == 3:
    with open('input.i3d', 'r') as file:
    
        # Read all lines into a list
        lines = file.readlines()
    
        # Extract the 28th line, where cpg option is specified 
        line = lines[27]  
    
        # Removing characters in front of the cpg option and the comment
        line = line.split('!')[0]
        cpg = line.split('=')[-1].strip()
                         
#!--- Reading of files section ---!

# Reading of friction coefficient history
M1 = np.loadtxt('cf_history.txt', skiprows=1, delimiter=',', dtype=np.float64)

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

    # Calculating the delta (in time units) between savings (skipping first saving)
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

    # Average
    mean_cf = np.mean(cfx[lower_index:])
    print("Mean cf value:", mean_cf)

#!--- Plot section, friction coefficient ---!

print("!--- Plotting of friction coefficient ---!")
print()

# Creating the folder for cf plot if it does not exist
if not os.path.exists("plots"):
    os.mkdir("plots")

# Subplots environment
fig, ax = plt.subplots(1, 1, figsize=(xinches,yinches), linewidth=tick_width, dpi=300)
   
# Friction coefficient
ax.scatter(time_unit, cfx, marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')

# Create a rectangle patch to show points we are excluding from average
rect = patches.Rectangle((0, 0), lower_tu, ylimsup, linewidth=0, edgecolor='none', facecolor='r', alpha=0.1)

# Channel flow only
if itype == 3:
    
    # Add the patch to the plot
    ax.add_patch(rect)
    
    # Horizontal line to show mean cf value
    ax.hlines(y=mean_cf, xmin=lower_tu, xmax=xlimsup, linewidth=lw, color=grey, linestyles='dashed', label=f'Mean value: {mean_cf:.3e}')
    
    # Creating the folder for cf average
    if not os.path.exists("data_post"):
        os.mkdir("data_post")
        
    




# Axes labels
ax.set_ylabel(r'$c_f$', fontsize=fla, labelpad=pad_axes_lab)

# x-axis label for Channel with CPG off
if itype == 3 and cpg == cpg_check:
    ax.set_xlabel(r'$t\frac{U_p}{h}$', fontsize=fla, labelpad=pad_axes_lab)
else:
    ax.set_xlabel(r'$t$', fontsize=fla, labelpad=pad_axes_lab)

# Axes limits
plt.xlim([xliminf, xlimsup])
plt.ylim([yliminf, ylimsup])

# Setting major and minor ticks parameters on both axes
ax.tick_params(axis='both', which='major', direction='in', length=lmajt, width=tick_width, top=True, right=True, pad=pad_numbers) 
ax.tick_params(axis='both', which='minor', direction='in', length=lmint, width=tick_width, top=True, right=True)

# Setting ticks
ax.tick_params(axis='both', labelcolor="k", labelsize=fla2)

# Saving the figure and show it
plt.savefig(f'plots/cf_vs_time{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)
plt.show()

#!--------------------------------------------------------------------------------------!








