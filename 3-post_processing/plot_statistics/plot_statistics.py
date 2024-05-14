#!---------------------------------------------------------!
#! With this script, we perform plotting of statistics     !
#! for TTBL simulations.                                   !
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

# Parameters
uwall = np.float64(1.0)    # wall velocity
re    = np.float64(500.0)  # Reynolds number
nu    = 1.0/re             # kinematic viscosity

#!--- Reading of files section ---!

# Reading of mean statistics
M = np.loadtxt('mean_stats600.0.txt', skiprows=1, delimiter=',', dtype=np.float64)

mean_u  = M[:,0]   
mean_v  = M[:,1]
var_u   = M[:,3]
var_v   = M[:,4]
mean_uv = M[:,12]

# Reading of vorticity components and mean gradient
M = np.loadtxt('vort_stats600.0.txt', skiprows=1, delimiter=',', dtype=np.float64)

vort_x = M[:,0]
vort_y = M[:,1]
vort_z = M[:,2]
mg     = M[:,3]

# Reading of grid points
y = np.loadtxt('yp.dat')

# Reading of the mean dissipation
M = np.loadtxt('diss_stats600.0.txt', skiprows=1, delimiter=',', dtype=np.float64)
eps = M[:]

#!--------------------------------!

#!--- Calculations ---!

# Shift due to the translating wall
mean_u = uwall - mean_u

# Shear quantities
sh_vel = np.sqrt(nu * mg[0])
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

# Viscous sub-layer & Von Karman law
y_plus_vsl = np.linspace(1, 15, 15)
u_plus_vsl = y_plus_vsl

k = 0.384
B = 4.173
y_plus_k = np.linspace(5, 180, 175)
u_plus_k = (1.0 / k) * np.log(y_plus_k) + B

# Kolmogorov time scale
tau_eta = (nu/eps)**0.5

#!--- Plot section, mean velocity profile ---!

# Labels and values for x-axis
labels = [r"$0.1$", r"$1$", r"$5$", r"$30$", r"$60$", r"$100$", r"$180$", r"$500$" ]
values = [0.1,    1.0, 5.0, 30.0, 60.0, 100.0, 180.0, 500.0]

lw = 1.5  # linewidth for plots
fig, ax = plt.subplots(1, 1, figsize=(14,10))

# Mean velocity profile plot
ax.scatter(y_plus, mean_u, color=blue, marker='o', linewidth=1.5, s=40, facecolors='none', edgecolors='C0')
ax.plot(y_plus_vsl, u_plus_vsl, color=grey, linestyle='--', linewidth=lw)
ax.plot(y_plus_k, u_plus_k, color=grey, linestyle='--', linewidth=lw)

# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=50, labelpad=20)
ax.set_ylabel(r'$U^+$', fontsize=50, labelpad=20)

# Legend
plt.legend(['Present', 'Viscous sublayer and log law (Kozul et al. (2016))'], loc='upper left', fontsize=18)

# Grid
plt.grid(True, linestyle='--')

plt.ylim([0, 20])

# Logarithmic x-axis
plt.semilogx()

# Setting x-ticks
ax.set_xticks(values, labels, color="k", size=20, rotation='horizontal')
ax.set_xticklabels(labels, fontsize=20, rotation=0, ha='center') 

# Setting y-ticks
ax.tick_params(axis='y', labelcolor="k", labelsize=20)

caption = 'Log law with constants: k = 0.384, B = 4.173 (Kozul et al. (2016))'
plt.text(0.11, 17.0, caption, horizontalalignment='left', verticalalignment='center', fontsize=16, fontweight='bold')

# Saving the figure
plt.savefig('umean.pdf', format='pdf', bbox_inches='tight')

plt.show()


#!--- Plot section, dissipation-related statistics ---!

# Labels and values for x-axis
labels = [r"$0.1$", r"$1$", r"$5$", r"$30$", r"$60$", r"$100$", r"$180$", r"$500$" ]
values = [0.1,    1.0, 5.0, 30.0, 60.0, 100.0, 180.0, 500.0]

lw = 1.5  # linewidth for plots
fig, ax = plt.subplots(1, 1, figsize=(14,10))

# Kolmogorov time-scale
ax.scatter(y_plus, tau_eta, color=blue, marker='o', linewidth=1.5, s=40, facecolors='none', edgecolors='C0')

# Axes labels
ax.set_xlabel(r'$y^+$',       fontsize=50, labelpad=20)
ax.set_ylabel(r'$\tau_\eta$', fontsize=50, labelpad=20)

# Grid
plt.grid(True, linestyle='--')

# Logarithmic x-axis
plt.semilogx()

# Setting x-ticks
ax.set_xticks(values, labels, color="k", size=20, rotation='horizontal')
ax.set_xticklabels(labels, fontsize=20, rotation=0, ha='center') 

# Setting y-ticks
ax.tick_params(axis='y', labelcolor="k", labelsize=20)

# Saving the figure
plt.savefig('kolmogorov_time_scale.pdf', format='pdf', bbox_inches='tight')

plt.show()





