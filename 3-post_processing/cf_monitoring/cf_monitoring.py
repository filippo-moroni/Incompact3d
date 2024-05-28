#!---------------------------------------------------------!
#! With this script, we perform plotting of friction       !
#! coefficient and its mean value calculation.             !
#!---------------------------------------------------------!

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

# Settings
np.seterr(divide='ignore', invalid='ignore')

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Computer Modern Sans serif",
})

plt.rcParams.update({'figure.autolayout': True})

# Parameters for plotting
lw          = 2.0             # linewidth for plots
markersize  = 80              # marker size for scatter plot
fla         = 80              # fontsize of labels of x and y axes (major labels, variables)
fla2        = 36              # fontsize of numbers of x and y axes 
xliminf     = 0.1             # x axis inferior limit
xalign      = xliminf*1.1     # value to adjust translation in x of captions
pad_numbers = 20              # pad of numbers on both axes
lmajt       = 30              # length of major ticks
lmint       = 15              # length of minor ticks
tick_width  = 1.5             # width of ticks and external box
y_location  = 0.75            # percentage of the y-axis limit for (automatic) positioning of captions

# Axes width
mpl.rcParams['axes.linewidth'] = tick_width

# Set some useful colors
grey = [0.5, 0.5, 0.5]

#!--------------------------------------------------------------------------------------!
              
#!--- Reading of files section ---!

print("!--- Plotting of friction coefficient ---!")

# Reading of friction coefficient history
M1 = np.loadtxt('cf_history.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
print()

# Extracting quantities from the full matrix
cfx       = M1[:,4] 
time_unit = M1[:,7] 

#!--------------------------------!

# Creating the folder for cf plot if it does not exist
if not os.path.exists("plots"):
    os.mkdir("plots")

#!--- Plot section, friction coefficient ---!

fig, ax = plt.subplots(1, 1, figsize=(14,10),linewidth=tick_width)
   
# Friction coefficient
ax.scatter(time_unit, cfx, marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')

plt.legend(['Fixed walls', 'Moving walls'], loc='upper left', fontsize=18)

# Axes labels
ax.set_xlabel(r'$T$', fontsize=fla, labelpad=20)
ax.set_ylabel(r'$c_f$', fontsize=fla, labelpad=20)

# Setting major and minor ticks on both axes
ax.tick_params(axis='both', which='major', direction='in', length=lmajt, width=tick_width, top=True, right=True, pad=pad_numbers) 
ax.tick_params(axis='both', which='minor', direction='in', length=lmint, width=tick_width, top=True, right=True)

# Setting y-ticks
ax.tick_params(axis='y', labelcolor="k", labelsize=fla2)

# Saving the figure and show it
plt.savefig('plots/cf_vs_time.pdf', format='pdf', bbox_inches='tight')
plt.show()

#!--------------------------------------------------------------------------------------!








