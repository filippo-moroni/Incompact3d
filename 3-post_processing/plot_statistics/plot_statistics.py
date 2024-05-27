#!---------------------------------------------------------!
#! With this script, we perform plotting of statistics     !
#! for TTBL and channel simulations:                       !
#!                                                         !
#! - mean statistics (mean[u], var[u], etc.)               !
#! - Kolmogorov time scale tau_eta                         !
#!---------------------------------------------------------!

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import matplotlib as mpl

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

# Axes width
mpl.rcParams['axes.linewidth'] = tick_width

# Set some useful colors
grey = [0.5, 0.5, 0.5]

# Parameter to switch between Lee & Moser reference or Cimarelli, 'Turbulence' lecture notes
iswitch = 1 # (0: Lee & Moser, 1: Cimarelli)

#!--------------------------------------------------------------------------------------!

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

#!--- Parameters & reference data ---!
if itype == 13:

    # TTBL
    uwall = np.float64(1.0)               # wall velocity
    re    = np.float64(500.0)             # Reynolds number
    nu    = 1.0/re                        # kinematic viscosity

elif itype == 3:
    
    # Channel
    re_cent = np.float64(4225.96)         # centerline Reynolds number of a laminar Poiseuille flow
    re_tau  = 0.123*(re_cent**0.875)      # corresponding estimated friction Reynolds number (Re_tau ~ 180)
    nu      = 1.0/re_cent                 # kinematic viscosity
               
    # Reading of Lee & Moser (2015) data
    M = np.loadtxt('data_lee_retau180.txt', skiprows=72, dtype=np.float64)
    y_plus_lm = M[:,1]
    mean_u_lm = M[:,2]
    
    M = np.loadtxt('data_lee_fluct_retau180.txt', skiprows=75, dtype=np.float64)
    var_u_lm   =   M[:,2]
    var_v_lm   =   M[:,3]
    mean_uv_lm = - M[:,5]
               
#!--- Reading of files section ---!

# Channel
if itype == 3:
    
    print("!--- Plotting of statistics for a channel ---!")

    # Reading of mean statistics
    M1 = np.loadtxt('data_post/mean_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
    # Reading of vorticity components and mean gradient
    M2 = np.loadtxt('data_post/vort_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)

# TTBL
elif itype == 13:

    print("!--- Plotting of statistics for a TTBL ---!")

    # Asking to the user the specific snapshot to show
    snap_numb = input("Enter the snapshot number to show (3 digits): ")
         
    # Reading of mean statistics
    M1 = np.loadtxt(f'data_post/mean_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
    # Reading of vorticity components and mean gradient
    M2 = np.loadtxt(f'data_post/vort_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)

print()

# Extracting quantities from the full matrices
mean_u  = M1[:,0]   
var_u   = M1[:,3]
var_v   = M1[:,4]
mean_uv = M1[:,12]

if itype == 3:
    mean_uv =  - mean_uv

vort_x = M2[:,0]
vort_y = M2[:,1]
vort_z = M2[:,2]
mg_tot = M2[:,3]
mg_x   = M2[:,4]
mg_z   = M2[:,5]

# Reading of grid points
y = np.loadtxt('yp.dat')

# Number of points in y direction
ny = len(y)

# Halve the points in y direction for a channel
if itype == 3:
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
sh_vel = np.sqrt(nu * np.abs(mg_x[0]))
delta_nu = nu / sh_vel
t_nu = nu / (sh_vel ** 2)

# Rescaling variables through wall units
y_plus   = y / delta_nu 
mean_u  /= sh_vel
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

#!--- Plot section, mean velocity profile, with selection dipending on the flow case ---!

# Mean velocity profile
fig, ax = plt.subplots(1, 1, figsize=(14,10),linewidth=tick_width)

# TTBL
if itype == 13:
    labels = [r"$0.1$", r"$1$", r"$5$", r"$10$", r"$30$", r"$60$", r"$100$", r"$180$", r"$500$"]
    values = [   0.1,      1.0,    5.0,    10.0,    30.0,    60.0,    100.0,    180.0,    500.0]
    xlimsup = 520.0
    ylimsup = 20.0
        
    # Mean velocity profile
    ax.scatter(y_plus, mean_u, marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')
    # Viscous sublayer and log law
    ax.plot(y_plus_vsl, u_plus_vsl, color=grey, linestyle='--', linewidth=lw)
    ax.plot(y_plus_k, u_plus_k, color=grey, linestyle='--', linewidth=lw)
    plt.legend(['Present', 'Viscous sublayer and log law'], loc='upper left', fontsize=18)
    
    # Caption
    caption = 'Log law with constants: k = 0.384, B = 4.173 (Kozul et al. (2016))'
    
# Channel    
elif itype == 3:
    labels = [r"$0.1$", r"$1$", r"$5$", r"$10$", r"$30$", r"$60$", r"$100$", r"$180$"]
    values = [   0.1,      1.0,    5.0,    10.0,    30.0,    60.0,    100.0,    180.0]
    xlimsup = 300.0
    ylimsup = 25.0
    
    # Mean velocity profile
    ax.scatter(y_plus[:ny], mean_u[:ny], marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')
    ax.scatter(y_plus_lm, mean_u_lm, marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C1')
    # Viscous sublayer and log law
    ax.plot(y_plus_vsl, u_plus_vsl, color=grey, linestyle='--', linewidth=lw)
    ax.plot(y_plus_k, u_plus_k, color=grey, linestyle='--', linewidth=lw)
    plt.legend(['Present', 'Lee and Moser (2015)', 'Viscous sublayer and log law'], loc='upper left', fontsize=18)
    
    # Caption
    if iswitch == 0:
        caption = 'Log law with constants: k = 0.384, B = 4.27 (Lee and Moser (2015))'
    elif iswitch == 1:
        caption  = 'Log law with constants: k = 0.37, B = 5.2 (Cimarelli, Turb. Lect. Notes)'
    caption2 = 'First points of Lee and Moser data not displayed' 
    
    # Plotting caption2
    plt.text(xalign, ylimsup - 6.5, caption2, fontsize=16, fontweight='bold', ha='left')

# Plotting caption
plt.text(xalign, ylimsup - 5.5, caption, fontsize=16, fontweight='bold', ha='left')

# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=fla, labelpad=20)
ax.set_ylabel(r'$U^+$', fontsize=fla, labelpad=20)

# Axes limits
plt.ylim([0, ylimsup])
plt.xlim([xliminf, xlimsup])

# Logarithmic x-axis and linear y-axis
ax.set_xscale('log')
ax.set_yscale('linear')

# Minor x-ticks based on log10
ax.xaxis.set_minor_locator(LogLocator(base=10,subs='all'))

# Setting x-ticks with values and labels
ax.set_xticks(values, labels, color="k", rotation='horizontal')
ax.set_xticklabels(labels, fontsize=fla2, rotation=0, ha='center')

# Setting major and minor ticks on both axes
ax.tick_params(axis='both', which='major', direction='in', length=lmajt, width=tick_width, top=True, right=True, pad=pad_numbers) 
ax.tick_params(axis='both', which='minor', direction='in', length=lmint, width=tick_width, top=True, right=True)

# Setting y-ticks
ax.tick_params(axis='y', labelcolor="k", labelsize=fla2)

# Saving the figure and show it
plt.savefig('umean.pdf', format='pdf', bbox_inches='tight')
plt.show()

#!--------------------------------------------------------------------------------------!

# <u'u'>
fig, ax = plt.subplots(1, 1, figsize=(14,10))

# TTBL
if itype == 13:
    labels = [r"$0.1$", r"$1$", r"$5$", r"$10$", r"$30$", r"$60$", r"$100$", r"$180$", r"$500$"]
    values = [   0.1,      1.0,    5.0,    10.0,    30.0,    60.0,    100.0,    180.0,    500.0]
    xlimsup = 520.0
    ylimsup = 8.0
    
    # <u'u'>
    ax.scatter(y_plus[:ny], var_u[:ny], marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')
    plt.legend(['Present'], loc='upper left', fontsize=18)
        
# Channel    
elif itype == 3:
    labels = [r"$0.1$", r"$1$", r"$5$", r"$10$", r"$30$", r"$60$", r"$100$", r"$180$"]
    values = [   0.1,      1.0,    5.0,    10.0,    30.0,    60.0,    100.0,    180.0]
    xlimsup = 300.0
    ylimsup = 8.0
    
    # <u'u'>
    ax.scatter(y_plus[:ny], var_u[:ny], marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')
    ax.scatter(y_plus_lm, var_u_lm, marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C1')
    plt.legend(['Present', 'Lee and Moser (2015)'], loc='upper left', fontsize=18)
    
    caption2 = 'First points of Lee and Moser data not displayed' 
    
    # Plotting caption2
    plt.text(xalign, ylimsup - 1.5, caption2, fontsize=16, fontweight='bold', ha='left')

# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=fla, labelpad=20)
ax.set_ylabel(r'$\langle u^{\prime 2} \rangle^+$', fontsize=fla, labelpad=20)

# Axes limits
plt.ylim([0, ylimsup])
plt.xlim([xliminf, xlimsup])

# Logarithmic x-axis and linear y-axis
ax.set_xscale('log')
ax.set_yscale('linear')

# Minor x-ticks based on log10
ax.xaxis.set_minor_locator(LogLocator(base=10,subs='all'))

# Setting x-ticks with values and labels
ax.set_xticks(values, labels, color="k", rotation='horizontal')
ax.set_xticklabels(labels, fontsize=fla2, rotation=0, ha='center')

# Setting major and minor ticks on both axes
ax.tick_params(axis='both', which='major', direction='in', length=lmajt, width=tick_width, top=True, right=True, pad=pad_numbers) 
ax.tick_params(axis='both', which='minor', direction='in', length=lmint, width=tick_width, top=True, right=True)

# Setting y-ticks
ax.tick_params(axis='y', labelcolor="k", labelsize=fla2)

# Saving the figure and show it
plt.savefig('uvar.pdf', format='pdf', bbox_inches='tight')
plt.show()

#!--------------------------------------------------------------------------------------!

# <v'v'>
fig, ax = plt.subplots(1, 1, figsize=(14,10))

# TTBL
if itype == 13:
    labels = [r"$0.1$", r"$1$", r"$5$", r"$10$", r"$30$", r"$60$", r"$100$", r"$180$", r"$500$"]
    values = [   0.1,      1.0,    5.0,    10.0,    30.0,    60.0,    100.0,    180.0,    500.0]
    xlimsup = 520.0
    ylimsup = 1.0
    
    # <v'v'>
    ax.scatter(y_plus[:ny], var_v[:ny], marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')
    plt.legend(['Present'], loc='upper left', fontsize=18)
        
# Channel    
elif itype == 3:
    labels = [r"$0.1$", r"$1$", r"$5$", r"$10$", r"$30$", r"$60$", r"$100$", r"$180$"]
    values = [   0.1,      1.0,    5.0,    10.0,    30.0,    60.0,    100.0,    180.0]
    xlimsup = 300.0
    ylimsup = 1.0
    
    # <v'v'>
    ax.scatter(y_plus[:ny], var_v[:ny], marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')
    ax.scatter(y_plus_lm, var_v_lm, marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C1')
    plt.legend(['Present', 'Lee and Moser (2015)'], loc='upper left', fontsize=18)
    
    caption2 = 'First points of Lee and Moser data not displayed' 
    
    # Plotting caption2
    plt.text(xalign, ylimsup - 0.18, caption2, fontsize=16, fontweight='bold', ha='left')

# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=fla, labelpad=20)
ax.set_ylabel(r'$\langle v^{\prime 2} \rangle^+$', fontsize=fla, labelpad=20)

# Axes limits
plt.ylim([0, ylimsup])
plt.xlim([xliminf, xlimsup])

# Logarithmic x-axis and linear y-axis
ax.set_xscale('log')
ax.set_yscale('linear')

# Minor x-ticks based on log10
ax.xaxis.set_minor_locator(LogLocator(base=10,subs='all'))

# Setting x-ticks with values and labels
ax.set_xticks(values, labels, color="k", rotation='horizontal')
ax.set_xticklabels(labels, fontsize=fla2, rotation=0, ha='center')

# Setting major and minor ticks on both axes
ax.tick_params(axis='both', which='major', direction='in', length=lmajt, width=tick_width, top=True, right=True, pad=pad_numbers) 
ax.tick_params(axis='both', which='minor', direction='in', length=lmint, width=tick_width, top=True, right=True)

# Setting y-ticks
ax.tick_params(axis='y', labelcolor="k", labelsize=fla2)

# Saving the figure and show it
plt.savefig('vvar.pdf', format='pdf', bbox_inches='tight')
plt.show()

#!--------------------------------------------------------------------------------------!

# <u'v'>
fig, ax = plt.subplots(1, 1, figsize=(14,10))

# TTBL
if itype == 13:
    labels = [r"$0.1$", r"$1$", r"$5$", r"$10$", r"$30$", r"$60$", r"$100$", r"$180$", r"$500$"]
    values = [   0.1,      1.0,    5.0,    10.0,    30.0,    60.0,    100.0,    180.0,    500.0]
    xlimsup = 520.0
    ylimsup = 1.0
    
    # <u'v'>
    ax.scatter(y_plus[:ny], mean_uv[:ny], marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')
    plt.legend(['Present'], loc='upper left', fontsize=18)
        
# Channel    
elif itype == 3:
    labels = [r"$0.1$", r"$1$", r"$5$", r"$10$", r"$30$", r"$60$", r"$100$", r"$180$"]
    values = [   0.1,      1.0,    5.0,    10.0,    30.0,    60.0,    100.0,    180.0]
    xlimsup = 300.0
    ylimsup = 1.0
    
    # <u'v'>
    ax.scatter(y_plus[:ny], mean_uv[:ny], marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')
    ax.scatter(y_plus_lm, mean_uv_lm, marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C1')
    plt.legend(['Present', 'Lee and Moser (2015)'], loc='upper left', fontsize=18)
    
    caption2 = 'First points of Lee and Moser data not displayed' 
    
    # Plotting caption2
    plt.text(xalign, ylimsup - 0.18, caption2, fontsize=16, fontweight='bold', ha='left')

# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=fla, labelpad=20)
ax.set_ylabel(r'$-\langle u^{\prime} v^{\prime}\rangle^+$', fontsize=fla, labelpad=20)

# Axes limits
plt.ylim([0, ylimsup])
plt.xlim([xliminf, xlimsup])

# Logarithmic x-axis and linear y-axis
ax.set_xscale('log')
ax.set_yscale('linear')

# Minor x-ticks based on log10
ax.xaxis.set_minor_locator(LogLocator(base=10,subs='all'))

# Setting x-ticks with values and labels
ax.set_xticks(values, labels, color="k", rotation='horizontal')
ax.set_xticklabels(labels, fontsize=fla2, rotation=0, ha='center')

# Setting major and minor ticks on both axes
ax.tick_params(axis='both', which='major', direction='in', length=lmajt, width=tick_width, top=True, right=True, pad=pad_numbers) 
ax.tick_params(axis='both', which='minor', direction='in', length=lmint, width=tick_width, top=True, right=True)

# Setting y-ticks
ax.tick_params(axis='y', labelcolor="k", labelsize=fla2)

# Saving the figure and show it
plt.savefig('uvmean.pdf', format='pdf', bbox_inches='tight')
plt.show()


#!--- Plot section, dissipation-related statistics ---!
# to do






