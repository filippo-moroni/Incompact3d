#!---------------------------------------------------------!
#! With this script, we perform plotting of statistics     !
#! for TTBL and channel simulations:                       !
#!                                                         !
#! - mean statistics (mean[u], var[u], etc.)               !
#! - Kolmogorov time scale tau_eta                         !
#!---------------------------------------------------------!

import numpy as np
import matplotlib.pyplot as plt

# Settings
np.seterr(divide='ignore', invalid='ignore')

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Computer Modern Sans serif",
})

plt.rcParams.update({'figure.autolayout': True})

# Set some useful colors
blue = [57 / 255, 106 / 255, 177 / 255]
grey = [0.5, 0.5, 0.5]

# Parameter to switch between Lee & Moser reference or Cimarelli, 'Turbulence' lecture notes
iswitch = 1 # (0: Lee & Moser, 1: Cimarelli)

# Read if we are plotting a channel or a TTBL
with open('input.i3d', 'r') as file:
    
    # Read all lines into a list
    lines = file.readlines()
    
    # Extract the 8th line, where itype is specified 
    eighth_line = lines[7]  
    
    # Removing characters in front of the itype value
    itype = eighth_line.split('=')[-1].strip()
    
    # Convert to integer
    itype = int(itype)

# Parameters
if itype == 13:
    uwall = np.float64(1.0)               # wall velocity
    re    = np.float64(500.0)             # Reynolds number
    nu    = 1.0/re                        # kinematic viscosity

elif itype == 3:
    
    re_cent = np.float64(4225.96)         # centerline Reynolds number of a laminar Poiseuille flow
    re_tau  = 0.123*(re_cent**0.875)      # corresponding estimated friction Reynolds number 
    nu      = 1.0/re_cent                 # kinematic viscosity
               
#!--- Reading of files section ---!

# Reading of mean statistics
M = np.loadtxt('data_post/mean_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)

mean_u  = M[:,0]   
mean_v  = M[:,1]
var_u   = M[:,3]
var_v   = M[:,4]
mean_uv = M[:,12]

# Reading of vorticity components and mean gradient
M = np.loadtxt('data_post/vort_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)

vort_x = M[:,0]
vort_y = M[:,1]
vort_z = M[:,2]
mg     = M[:,3]

# Reading of grid points
y = np.loadtxt('yp.dat')

# Number of points in y direction
ny = len(y)
ny = (ny - 1) // 2 + 1

# Reading of the mean dissipation
#M = np.loadtxt('data_post/diss_stats-030.txt', skiprows=1, delimiter=',', dtype=np.float64)
#eps = M[:]

#!--------------------------------!

#!--- Calculations ---!

# Valid only for TTBLs
if itype == 13:

    # Shift due to the translating wall
    mean_u = uwall - mean_u
           
# Shear quantities
sh_vel = np.sqrt(nu * mg[0])
delta_nu = nu / sh_vel
t_nu = nu / (sh_vel ** 2)

# Rescaling variables through wall units
y_plus   = y / delta_nu 
mean_u   = mean_u / sh_vel
var_u   /= sh_vel ** 2
var_v   /= sh_vel ** 2
mean_uv /= sh_vel ** 2

vort_x *= t_nu
vort_y *= t_nu
vort_z *= t_nu

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

    if iswitch == 0:
    
        # Lee and Moser (2015)
        k = 0.384
        B = 4.27
    
    elif iswitch == 1:
        
        # Cimarelli ('Turbulence' lecture notes)
        k = 0.37
        B = 5.2

# Von Karman law
y_plus_k = np.linspace(5, 180, 175)
u_plus_k = (1.0 / k) * np.log(y_plus_k) + B

#!-------------------------------!

# Kolmogorov time scale
#tau_eta = np.sqrt(nu/eps)

#!--- Plot section, mean velocity profile ---!

# Labels and values for x-axis
labels = [r"$0.1$", r"$1$", r"$5$", r"$30$", r"$60$", r"$100$", r"$180$", r"$500$" ]
values = [0.1,    1.0, 5.0, 30.0, 60.0, 100.0, 180.0, 500.0]

lw = 1.5  # linewidth for plots
fig, ax = plt.subplots(1, 1, figsize=(14,10))

# Mean velocity profile plot
if itype == 13:
    ax.scatter(y_plus, mean_u, color=blue, marker='o', linewidth=1.5, s=40, facecolors='none', edgecolors='C0')
elif itype == 3:
    ax.scatter(y_plus[:ny], mean_u[:ny], color=blue, marker='o', linewidth=1.5, s=40, facecolors='none', edgecolors='C0')

ax.plot(y_plus_vsl, u_plus_vsl, color=grey, linestyle='--', linewidth=lw)
ax.plot(y_plus_k, u_plus_k, color=grey, linestyle='--', linewidth=lw)

# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=50, labelpad=20)
ax.set_ylabel(r'$U^+$', fontsize=50, labelpad=20)

# Legend
plt.legend(['Present', 'Viscous sublayer and log law'], loc='upper left', fontsize=18)

# Grid
plt.grid(True, linestyle='--')

# y-axis limit
ylim = 22.0
plt.ylim([0, ylim])

# Logarithmic x-axis
ax.semilogx()

# Setting x-ticks
ax.set_xticks(values, labels, color="k", size=20, rotation='horizontal')
ax.set_xticklabels(labels, fontsize=20, rotation=0, ha='center') 

# Setting y-ticks
ax.tick_params(axis='y', labelcolor="k", labelsize=20)

# Differencing caption based on flow case
if itype == 13:
    # TTBL
    caption = 'Log law with constants: k = 0.384, B = 4.173 (Kozul et al. (2016))'
elif itype == 3:
    # Channel
    if iswitch == 0:
        caption = 'Log law with constants: k = 0.384, B = 4.27 (Lee and Moser (2015))'
    elif iswitch == 1:
        caption = 'Log law with constants: k = 0.37, B = 5.2 (Cimarelli, Turb. Lect. Notes)'

plt.text(0.12, ylim - 3.5, caption, horizontalalignment='left', verticalalignment='center', fontsize=16, fontweight='bold')

# Saving the figure
plt.savefig('umean.pdf', format='pdf', bbox_inches='tight')

plt.show()


#!--- Plot section, dissipation-related statistics ---!

# Index up to which you want to plot
#index = 50

# Labels and values for x-axis
#labels = [r"$0.1$", r"$1$", r"$5$", r"$30$", r"$60$", r"$100$", r"$180$", r"$500$" ]
#values = [0.1,    1.0, 5.0, 30.0, 60.0, 100.0, 180.0, 500.0]

#lw = 1.5  # linewidth for plots
#fig, ax = plt.subplots(1, 1, figsize=(14,10))

# Kolmogorov time-scale
#ax.scatter(y_plus[:index], tau_eta[:index], color=blue, marker='o', linewidth=1.5, s=40, facecolors='none', edgecolors='C0')

# Axes labels
#ax.set_xlabel(r'$y^+$',       fontsize=50, labelpad=20)
#ax.set_ylabel(r'$\tau_\eta$', fontsize=50, labelpad=20)

# Grid
#plt.grid(True, linestyle='--')

# Logarithmic x-axis
#plt.semilogx()

# Setting x-ticks
#ax.set_xticks(values, labels, color="k", size=20, rotation='horizontal')
#ax.set_xticklabels(labels, fontsize=20, rotation=0, ha='center') 

# Setting y-ticks
#ax.tick_params(axis='y', labelcolor="k", labelsize=20)

# Saving the figure
#plt.savefig('kolmogorov_time_scale.pdf', format='pdf', bbox_inches='tight')

#plt.show()





