#!---------------------------------------------------------!
#! With this script, we perform plotting of statistics     !
#! for channel flow simulations, fixed walls compared to   !
#! oscillating walls.                                      !
#!                                                         !
#! - mean statistics (mean[u], var[u], etc.)               !
#!---------------------------------------------------------!

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import matplotlib as mpl
import os

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
xliminf      = 0.1             # x axis inferior limit (y+)

# Page settings (A4 paper format: 8.3 x 11.7 inches)
xinches      = 2.6             # size in inches in x direction of the image
yinches      = 2.2             # size in inches in y direction of the image

# Axes width
mpl.rcParams['axes.linewidth'] = tick_width

# Set some useful colors
grey = [0.5, 0.5, 0.5]

# Parameter to switch between Lee & Moser reference or Cimarelli, 'Turbulence' lecture notes
iswitch = 1 # (0: Lee & Moser, 1: Cimarelli)

#!--------------------------------------------------------------------------------------!

# Read the name of the flowcase
with open('post.prm', 'r') as file:
    
    # Read all lines into a list
    lines = file.readlines()
    
    # Extract flowcase name
    add_string = lines[3] 
    add_string = add_string.split('!')[0]
    add_string = add_string.rstrip()

# Read input file
with open('input.i3d', 'r') as file:
    
    # Read all lines into a list
    lines = file.readlines()
    
    # Extract itype, nx, nz, Lx, Ly, Lz, Re 
    itype = lines[7]  
    nx    = lines[14]
    nz    = lines[16]
    Lx    = lines[21]
    Ly    = lines[22]
    Lz    = lines[23]
    re    = lines[26]
    
    # Removing characters in front of the extracted strings and the comments
    itype = itype.split('=')[-1].strip()
    
    nx    = nx.split('!')[0]
    nx    = nx.split('=')[-1].strip()
    
    nz    = nz.split('!')[0]
    nz    = nz.split('=')[-1].strip()
    
    Lx    = Lx.split('!')[0]
    Lx    = Lx.split('=')[-1].strip()
    
    Ly    = Ly.split('!')[0]
    Ly    = Ly.split('=')[-1].strip()
    
    Lz    = Lz.split('!')[0]
    Lz    = Lz.split('=')[-1].strip()
    
    re    = re.split('!')[0]
    re    = re.split('=')[-1].strip()
    
    # Convert to integer
    itype = int(itype)
    nx    = int(nx)
    nz    = int(nz)
    Lx    = np.float64(Lx)
    Ly    = np.float64(Ly)
    Lz    = np.float64(Lz)
    re    = np.float64(re)
    
#!--- Parameters & reference data ---!

# Channel (valid only for CFR at the moment)
re_cent = np.float64(re)             # centerline Reynolds number of a laminar Poiseuille flow
re_tau  = 0.123*(re_cent**0.875)     # corresponding estimated friction Reynolds number 
nu      = 1.0/re_cent                # kinematic viscosity
               
# Reading of Lee & Moser (2015) data
M = np.loadtxt('reference_data/data_lee_retau180.txt', skiprows=72, dtype=np.float64)
y_plus_lm = M[:,1]
mean_u_lm = M[:,2]
    
M = np.loadtxt('reference_data/data_lee_fluct_retau180.txt', skiprows=75, dtype=np.float64)
var_u_lm   =   M[:,2]
var_v_lm   =   M[:,3]
mean_uv_lm = - M[:,5]
               
#!--- Reading of files section ---!
print()

print("!--- Plotting of statistics for a channel ---!")

# Reading of mean statistics
M1 = np.loadtxt('data_post/mean_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
# Reading of vorticity components and mean gradient
M2 = np.loadtxt('data_post/vort_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)

print()

# Extracting quantities from the full matrices
mean_u  = M1[:,0]   
var_u   = M1[:,3]
var_v   = M1[:,4]
mean_uv = M1[:,12]

# Change sign for Reynolds stresses
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

#!--------------------------------!

#!--- Calculations ---!

# Mesh spacings
delta_x = Lx / nx
delta_z = Lz / nz
           
# Shear quantities
sh_vel = np.sqrt(nu * np.abs(mg_x[0]))
delta_nu = nu / sh_vel
t_nu = nu / (sh_vel ** 2)

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

vort_x *= t_nu
vort_y *= t_nu
vort_z *= t_nu

#!--- Reference mean profiles ---!

# Viscous sub-layer
y_plus_vsl = np.linspace(1, 15, 15)
u_plus_vsl = y_plus_vsl

# Log law constants 
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

# Creating the folder for plots if it does not exist
os.makedirs('plots', mode=0o777, exist_ok=True)
    
#!--- Plot section, mean velocity profile ---!

# Mean velocity profile
fig, ax = plt.subplots(1, 1, figsize=(xinches,yinches), linewidth=tick_width, dpi=300)

xlimsup = 300.0
ylimsup = 25.0
    
# Mean velocity profile
ax.scatter(y_plus[:ny], mean_u[:ny], marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')
ax.plot(y_plus_lm, mean_u_lm, color='C1', linestyle='-', linewidth=lw)
    
# Viscous sublayer and log law
ax.plot(y_plus_vsl, u_plus_vsl, color=grey, linestyle='--', linewidth=lw)
ax.plot(y_plus_k, u_plus_k, color=grey, linestyle='--', linewidth=lw)
    
# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=fla, labelpad=pad_axes_lab)
ax.set_ylabel(r'$U^+$', fontsize=fla, labelpad=pad_axes_lab)

# Axes limits
plt.ylim([0, ylimsup])
plt.xlim([xliminf, xlimsup])

# Logarithmic x-axis and linear y-axis
ax.set_xscale('log')
ax.set_yscale('linear')

# Minor x-ticks based on log10
ax.xaxis.set_minor_locator(LogLocator(base=10,subs='all'))
    
# Setting major and minor ticks on both axes
ax.tick_params(axis='both', which='major', direction='in', length=lmajt, width=tick_width, pad=pad_numbers, labelsize=fla2, labelcolor='k') 
ax.tick_params(axis='both', which='minor', direction='in', length=lmint, width=tick_width)

# Setting x-ticks labels
plt.xticks(ha='left')

# Saving the figure
plt.savefig(f'plots/umean_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)
    
# Show the figure
plt.show()

#!--------------------------------------------------------------------------------------!

# <u'u'>
fig, ax = plt.subplots(1, 1, figsize=(xinches,yinches), linewidth=tick_width, dpi=300)

xlimsup = 300.0
ylimsup = 8.0
    
# Plotting
ax.scatter(y_plus[:ny], var_u[:ny], marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')
ax.plot(y_plus_lm, var_u_lm, color='C1', linestyle='-', linewidth=lw)
    
# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=fla, labelpad=pad_axes_lab)
ax.set_ylabel(r'$\langle u^{\prime 2} \rangle^+$', fontsize=fla, labelpad=pad_axes_lab)

# Axes limits
plt.ylim([0, ylimsup])
plt.xlim([xliminf, xlimsup])

# Logarithmic x-axis and linear y-axis
ax.set_xscale('log')
ax.set_yscale('linear')

# Minor x-ticks based on log10
ax.xaxis.set_minor_locator(LogLocator(base=10,subs='all'))
    
# Setting major and minor ticks on both axes
ax.tick_params(axis='both', which='major', direction='in', length=lmajt, width=tick_width, pad=pad_numbers, labelsize=fla2, labelcolor='k') 
ax.tick_params(axis='both', which='minor', direction='in', length=lmint, width=tick_width)

# Setting x-ticks labels
plt.xticks(ha='left')

# Saving the figure
plt.savefig(f'plots/uvar_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)

# Show the figure
plt.show()

#!--------------------------------------------------------------------------------------!

# <v'v'>
fig, ax = plt.subplots(1, 1, figsize=(xinches,yinches), linewidth=tick_width, dpi=300)

xlimsup = 300.0
ylimsup = 0.8
    
# Plotting
ax.scatter(y_plus[:ny], var_v[:ny], marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')
ax.plot(y_plus_lm, var_v_lm, color='C1', linestyle='-', linewidth=lw)
    
# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=fla, labelpad=pad_axes_lab)
ax.set_ylabel(r'$\langle v^{\prime 2} \rangle^+$', fontsize=fla, labelpad=pad_axes_lab)

# Axes limits
plt.ylim([0, ylimsup])
plt.xlim([xliminf, xlimsup])

# Logarithmic x-axis and linear y-axis
ax.set_xscale('log')
ax.set_yscale('linear')

# Minor x-ticks based on log10
ax.xaxis.set_minor_locator(LogLocator(base=10,subs='all'))
    
# Setting major and minor ticks on both axes
ax.tick_params(axis='both', which='major', direction='in', length=lmajt, width=tick_width, pad=pad_numbers, labelsize=fla2, labelcolor='k') 
ax.tick_params(axis='both', which='minor', direction='in', length=lmint, width=tick_width)

# Setting x-ticks labels
plt.xticks(ha='left')

# Saving the figure
plt.savefig(f'plots/vvar_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)

# Show the figure
plt.show()

#!--------------------------------------------------------------------------------------!

# <u'v'>
fig, ax = plt.subplots(1, 1, figsize=(xinches,yinches), linewidth=tick_width, dpi=300)

xlimsup = 300.0
ylimsup = 0.8

# Plotting    
ax.scatter(y_plus[:ny], mean_uv[:ny], marker='o', linewidth=lw, s=markersize, facecolors='none', edgecolors='C0')
ax.plot(y_plus_lm, mean_uv_lm, color='C1', linestyle='-', linewidth=lw)
       
# y-axis label
ax.set_ylabel(r'$-\langle u^{\prime} v^{\prime}\rangle^+$', fontsize=fla, labelpad=pad_axes_lab)

# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=fla, labelpad=pad_axes_lab)

# Axes limits
plt.ylim([0, ylimsup])
plt.xlim([xliminf, xlimsup])

# Logarithmic x-axis and linear y-axis
ax.set_xscale('log')
ax.set_yscale('linear')

# Minor x-ticks based on log10
ax.xaxis.set_minor_locator(LogLocator(base=10,subs='all'))
    
# Setting major and minor ticks on both axes
ax.tick_params(axis='both', which='major', direction='in', length=lmajt, width=tick_width, pad=pad_numbers, labelsize=fla2, labelcolor='k') 
ax.tick_params(axis='both', which='minor', direction='in', length=lmint, width=tick_width)

# Setting x-ticks labels
plt.xticks(ha='left')

# Saving the figure
plt.savefig(f'plots/uvmean_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)

# Show the figure
plt.show()







