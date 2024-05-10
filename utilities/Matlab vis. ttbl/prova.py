import numpy as np
import matplotlib.pyplot as plt

# Input from the user
Re = 500.0                # Reynolds number (1/nu)
nu = 1.0/Re               # kinematic viscosity

ny = 649                  # number of points in y-direction
Ly = 65.1                 # total height of the domain

# Reading of file and variables
M1 = pd.read_csv('mean_stats430.0.txt', skiprows=1, delim_whitespace=True)
M3 = pd.read_csv('vort_stats430.0.txt', skiprows=1, delim_whitespace=True)
G = pd.read_csv('yp.dat', header=None)

# Default code variables
mean_u = M1.iloc[:, 0].values        # mean of u default code
mean_v = M1.iloc[:, 1].values        # mean of v default code
var_u = M1.iloc[:, 3].values         # variance of u
var_v = M1.iloc[:, 4].values         # variance of v
mean_uv = M1.iloc[:, 12].values      # <u'v'>

uwall = 1.0
mean_u = uwall - mean_u

# Vorticity
vort_x = M3.iloc[:, 0].values
vort_y = M3.iloc[:, 1].values
vort_z = M3.iloc[:, 2].values

y = G.iloc[:, 0].values

# Calculations for default code
mean_gradient = mean_u[1] / y[1]     # partial U / partial y
sh_vel = np.sqrt(nu * mean_gradient)
delta_nu = nu / sh_vel
t_nu = nu / (sh_vel ** 2)

# Rescaling variables through wall units
y_plus = y / delta_nu
mean_u = mean_u / sh_vel
var_u = var_u / (sh_vel ** 2)
var_v = var_v / (sh_vel ** 2)
mean_uv = mean_uv / (sh_vel ** 2)

vort_x = vort_x * t_nu
vort_y = vort_y * t_nu
vort_z = vort_z * t_nu

# Viscous sub-layer & Von Karman law
y_plus_vsl = np.linspace(1, 15, 15)
u_plus_vsl = y_plus_vsl

k = 0.384
B = 4.173
y_plus_k = np.linspace(5, 180, 175)
u_plus_k = (1 / k) * np.log(y_plus_k) + B

# Mean velocity profile plot
plt.figure(figsize=(22, 12))

plt.scatter(y_plus, mean_u, marker='o', edgecolors='blue', linewidth=1.5)
plt.plot(y_plus_vsl, u_plus_vsl, color='grey', linestyle='--', linewidth=1.5)
plt.plot(y_plus_k, u_plus_k, color='grey', linestyle='--', linewidth=1.5)

plt.legend(['Present', 'Viscous sublayer and log law (Kozul et al. (2016))'], loc='upper left', fontsize=18)

plt.grid(True, which='both')
plt.minorticks_on()

plt.xlim([0, 180])
plt.xscale('log')
plt.xticks([0, 5, 30, 60, 100, 180])
plt.xlabel('$y^+$', fontsize=50)

yaxis_lim = 20
plt.ylim([0, yaxis_lim])
plt.ylabel('$U^+$', fontsize=50)

plt.text(0.32, -1, 'Log law with constants: k = 0.384, B = 4.173 (Kozul et al. (2016))', 
         horizontalalignment='left', verticalalignment='middle', fontsize=16, fontweight='bold')

plt.show()

# Vorticity plot
'''
plt.figure(figsize=(22, 12))

plt.scatter(y_plus, vort_z, marker='o', edgecolors='blue', linewidth=1.5)

plt.legend(['Present'], loc='upper left', fontsize=12)

plt.grid(True, which='both')
plt.minorticks_on()

plt.xlim([0, y_plus[ny]])
plt.xscale('log')
plt.xticks([0, 1, 5, 30, 60, 100, 180, 300, y_plus[ny]])
plt.xlabel('$y^+$', fontsize=40)

yaxis_lim = 20
plt.ylabel('$\langle \omega_z \rangle$', fontsize=40)

plt.show()
'''

