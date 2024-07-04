#!---------------------------------------------------------!
#! With this script, we perform plotting of statistics     !
#! for TTBLs and channel flow simulations:                 !
#!                                                         !
#! - mean statistics (mean[u], var[u], etc.);              !
#! - mean total dissipation (to be done);                  !
#! - correlation coefficients for spanwise correlations.   !
#!                                                         !
#! Calculated and stored:                                  !
#!                                                         !
#! - non-dimensional grid spacings and domain dimensions.  !
#!---------------------------------------------------------!

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

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
from read_input_files import read_input_files

#!--------------------------------------------------------------------------------------!

# Create folders to store later results (e.g. grid spacings and time scales files, plots)
os.makedirs('data_post', mode=0o777, exist_ok=True)
os.makedirs('plots',     mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
itype, nx, nz, Lx, Ly, Lz, re, iswitch_wo, add_string = read_input_files('input.i3d','post.prm')

# Now you can use the variables as needed
print(f"itype: {itype}")
print(f"nx: {nx}")
print(f"nz: {nz}")
print(f"Lx: {Lx}")
print(f"Ly: {Ly}")
print(f"Lz: {Lz}")
print(f"Re: {re}")
print(f"iswitch_wo: {iswitch_wo}")

exit()


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
    
    # Extract itype, nx, nz, Lx, Ly, Lz, Re, iswitch_wo 
    itype      = lines[7]  
    nx         = lines[14]
    nz         = lines[16]
    Lx         = lines[21]
    Ly         = lines[22]
    Lz         = lines[23]
    re         = lines[26]
    iswitch_wo = lines[91]
    
    # Removing characters in front of the extracted strings and the comments:
    # 1) split: the string is split when the specified character is encountered; 
    # 2) we select the portion of string with index inside square brackets;
    # 3) strip: removes leading or trailing whitespaces from the string. 
    
    itype = itype.split('=')[-1].strip()
    
    nx         = nx.split('!')[0]
    nx         = nx.split('=')[-1].strip()
    
    nz         = nz.split('!')[0]
    nz         = nz.split('=')[-1].strip()
    
    Lx         = Lx.split('!')[0]
    Lx         = Lx.split('=')[-1].strip()
    
    Ly         = Ly.split('!')[0]
    Ly         = Ly.split('=')[-1].strip()
    
    Lz         = Lz.split('!')[0]
    Lz         = Lz.split('=')[-1].strip()
    
    re         = re.split('!')[0]
    re         = re.split('=')[-1].strip()
    
    iswitch_wo = iswitch_wo.split('!')[0]
    iswitch_wo = iswitch_wo.split('=')[-1].strip()
    
    # Convert to needed variable type (integer, float, etc.)
    itype      = int(itype)
    nx         = int(nx)
    nz         = int(nz)
    Lx         = np.float64(Lx)
    Ly         = np.float64(Ly)
    Lz         = np.float64(Lz)
    re         = np.float64(re)
    iswitch_wo = int(iswitch_wo)
    
#!--- Parameters & reference data ---!
if itype == 13:

    # TTBL
    uwall = np.float64(1.0)              # wall velocity
    re    = np.float64(re)               # Reynolds number
    nu    = 1.0/re                       # kinematic viscosity

elif itype == 3:
    
    # Channel (valid only for CFR at the moment)
    re_cent = np.float64(re)             # centerline Reynolds number of a laminar Poiseuille flow
    re_tau  = 0.123*(re_cent**0.875)     # corresponding estimated friction Reynolds number 
    nu      = 1.0/re_cent                # kinematic viscosity
               
    # Reading of Lee & Moser (2015) data
    M = np.loadtxt('reference_data/lee&moser2015/mean_stats_lee&moser2015.txt', skiprows=72, dtype=np.float64)
    y_plus_lm = M[:,1]
    mean_u_lm = M[:,2]
    
    M = np.loadtxt('reference_data/lee&moser2015/var_stats_lee&moser2015.txt', skiprows=75, dtype=np.float64)
    var_u_lm   =   M[:,2]
    var_v_lm   =   M[:,3]
    mean_uv_lm = - M[:,5]
    
    # Velocity auto-correlations, Kim et al. (1987) data, y+ = 10.52
    M = np.loadtxt('reference_data/kim1987/cuuz_kim1987.txt', skiprows=7, delimiter=',', dtype=np.float64)
    rz_plus_cuuz_kim = M[:,0]
    cuuz_kim         = M[:,1]
    
    M = np.loadtxt('reference_data/kim1987/cvvz_kim1987.txt', skiprows=7, delimiter=',', dtype=np.float64)
    rz_plus_cvvz_kim = M[:,0]
    cvvz_kim         = M[:,1]
    
    M = np.loadtxt('reference_data/kim1987/cwwz_kim1987.txt', skiprows=7, delimiter=',', dtype=np.float64)
    rz_plus_cwwz_kim = M[:,0]
    cwwz_kim         = M[:,1]
        
    # Reading of wall-oscillations data (A^+ = 12, T^+ = 100) 
    if iswitch_wo == 1:
    
        # Mean velocity profile (Touber & Leschziner (2012))
        M = np.loadtxt('reference_data/touber2012/umean_touber2012.txt', skiprows=8, delimiter=',', dtype=np.float64)
        y_plus_touber = M[:,0]
        mean_u_touber = M[:,1]
        
        # Mean velocity profile, Reynolds stress and RMSs (Yao et al. (2019))
        M = np.loadtxt('reference_data/yao2019/umean_yao2019.txt', skiprows=8, delimiter=',', dtype=np.float64)
        y_plus_umean_yao = M[:,0]
        mean_u_yao       = M[:,1]
        
        M = np.loadtxt('reference_data/yao2019/uvar_yao2019.txt', skiprows=8, delimiter=',', dtype=np.float64)
        y_plus_uvar_yao = M[:,0]
        var_u_yao       = M[:,1]
        
        M = np.loadtxt('reference_data/yao2019/vvar_yao2019.txt', skiprows=8, delimiter=',', dtype=np.float64)
        y_plus_vvar_yao = M[:,0]
        var_v_yao       = M[:,1]
        
        M = np.loadtxt('reference_data/yao2019/uvmean_yao2019.txt', skiprows=8, delimiter=',', dtype=np.float64)
        y_plus_uvmean_yao = M[:,0]
        mean_uv_yao       = M[:,1]
        
        # Rescale Yao et al. (2019) data from uncontrolled to controlled shear velocity (0: uncontrolled, c: controlled)
        cf_0_yao = 0.00790
        cf_c_yao = 0.00511
        
        sh_vel_0_yao = (2.0/3.0)*np.sqrt(cf_0_yao / 2.0)
        sh_vel_c_yao = (2.0/3.0)*np.sqrt(cf_c_yao / 2.0)
        
        # Mean velocity is already rescaled by the actual shear velocity in Yao's data
        
        # Rescale RMSs and obtain variances
        y_plus_uvar_yao   = y_plus_uvar_yao   * (sh_vel_c_yao / sh_vel_0_yao)
        y_plus_vvar_yao   = y_plus_vvar_yao   * (sh_vel_c_yao / sh_vel_0_yao)
        y_plus_uvmean_yao = y_plus_uvmean_yao * (sh_vel_c_yao / sh_vel_0_yao)
        
        var_u_yao   = (var_u_yao   *  sh_vel_0_yao / sh_vel_c_yao)**2
        var_v_yao   = (var_v_yao   *  sh_vel_0_yao / sh_vel_c_yao)**2
        mean_uv_yao = mean_uv_yao  * (sh_vel_0_yao / sh_vel_c_yao)**2
        
                                                                              
#!--- Reading of files section ---!
print()

# Channel
if itype == 3:
    
    print("!--- Plotting of statistics for a channel ---!")

    # Reading of mean statistics
    M1 = np.loadtxt('data_post/mean_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
    # Reading of vorticity components and mean gradient
    M2 = np.loadtxt('data_post/vort_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
    # Reading of the mean total dissipation
    eps = np.loadtxt('data_post/diss_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
    # Reading of correlations
    Ruuz = np.loadtxt('data_post/Ruuz.txt', skiprows=0, delimiter=None, dtype=np.float64)
    Rvvz = np.loadtxt('data_post/Rvvz.txt', skiprows=0, delimiter=None, dtype=np.float64)
    Rwwz = np.loadtxt('data_post/Rwwz.txt', skiprows=0, delimiter=None, dtype=np.float64)
        
# TTBL
elif itype == 13:

    print("!--- Plotting of statistics for a TTBL ---!")

    # Asking to the user the specific snapshot to show
    snap_numb = input("Enter the snapshot number to show (4 digits): ")
         
    # Reading of mean statistics
    M1 = np.loadtxt(f'data_post/mean_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
    # Reading of vorticity components and mean gradient
    M2 = np.loadtxt(f'data_post/vort_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
    # Reading of the mean total dissipation
    eps = np.loadtxt(f'data_post/diss_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
        
    # Reading of correlations
    Ruuz = np.loadtxt(f'data_post/Ruuz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)
    Rvvz = np.loadtxt(f'data_post/Rvvz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)
    Rwwz = np.loadtxt(f'data_post/Rwwz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)

print()

# Extracting quantities from the full matrices
mean_u  = M1[:,0]
mean_w  = M1[:,2]   
var_u   = M1[:,3]
var_v   = M1[:,4]
mean_uv = M1[:,12]

# Valid only for Channel
if itype == 3:
    
    # Change sign for Reynolds stresses
    mean_uv =  - mean_uv

# Valid only for TTBLs
elif itype == 13:

    # Shift due to the translating wall
    mean_u = uwall - mean_u

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

# Select the height at which correlations are plotted
y_plus_in = np.float64(input("Enter y+ value for correlations plotting: "))
print()

#!--------------------------------------------------------------------------------------!

#!--- Calculations ---!

# Mesh spacings
delta_x = Lx / nx
delta_z = Lz / nz
           
# Shear quantities
sh_vel   = np.sqrt(nu * np.abs(mg_x[0]))  # shear velocity (based on streamwise mean gradient)  
delta_nu = nu / sh_vel                    # viscous length
t_nu     = nu / (sh_vel ** 2)             # viscous time

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

# Print viscous time unit
print("Viscous time unit, t_nu = ", t_nu)
print()

# y+ at the centerline or at the BL edge
if itype == 3:
    delta_yd_plus = y_plus[ny-1] - y_plus[ny-2]
elif itype == 13:
    # Calculate the index at which the BL thickness delta99 is
    j = 0
    while mean_u[j] > mean_u[0]*0.01: j = j + 1
    delta_yd_plus = y_plus[j] - y_plus[j-1] 

#!--------------------------------------------------------------------------------------!
    
#!--- Writing to file the non-dimensional grid spacings and domain dimensions ---!
           
# Create the file and write  
with open('data_post/grid_spacings_post.txt', 'w') as f:
    f.write(f"{'delta_x^+':<{pp.c_w}}, "  +
            f"{'delta_y1^+':<{pp.c_w}}, " +
            f"{'delta_z^+':<{pp.c_w}}, "  +
            f"{'Lx^+':<{pp.c_w}}, "       +
            f"{'Ly^+/2':<{pp.c_w}}, "     +
            f"{'Lz^+':<{pp.c_w}}, "       +
            f"{'delta_yd^+':<{pp.c_w}}\n" )

    f.write(f"{delta_x_plus:{pp.fs}}, "   +
            f"{y_plus[1]:{pp.fs}}, "      +
            f"{delta_z_plus:{pp.fs}}, "   +
            f"{Lx_plus:{pp.fs}}, "        +
            f"{Ly_plus/2:{pp.fs}}, "      +
            f"{Lz_plus:{pp.fs}}, "        +
            f"{delta_yd_plus:{pp.fs}}\n"  ) 

#!--------------------------------------------------------------------------------------!

# Find the maximum of mean total dissipation
eps_max = max(eps)

# Minimum Kolmogorov time scale and print it
tau_eta = np.sqrt(nu/eps_max)
print("Minimum Kolmogorov time scale, tau_eta = ", tau_eta)
print()

#!--- Calculations for correlations ---!

# Search for the index corresponding to the target y+ for correlations
c = 0 
for j in range(0, ny-1, 1):   
    if y_plus[j] < y_plus_in: c = c + 1

# Print the actual y+ value selected
print("Actual y+ value selected = ", y_plus[c])
print()

# Take the Rii value at rz = 0 and rescale to obtain correlation coefficients
temp = Ruuz[c,0]
Ruuz = Ruuz / temp

temp = Rvvz[c,0]
Rvvz = Rvvz / temp

temp = Rwwz[c,0]
Rwwz = Rwwz / temp

# Halve the number of points in z-dir. to avoid periodicity effects
nz = nz // 2
Lz = Lz / 2.0

# Create the separation variable array
rz = np.linspace(0, Lz, nz)

# Calculate the index at which the correlation coefficient goes to zero
k = 0
while Ruuz[c,k] > 0.0: k = k + 1

rz0    = rz[0]    # First element of rz vector (rz = 0)
rzstar = rz[k-1]  # Element of rz vector at which Cii(rz) goes to zero

# Rescale separation variable by viscous unit
rz = rz / delta_nu
    
#!--- Calculate the integral length scale lambda z ---!

# Interpolation at the 6th order of accuracy with a spline of 5th order
spl = InterpolatedUnivariateSpline(rz[:k], Ruuz[c,:k], k=5)
lambda_z = spl.integral(rz0, rzstar)

# Rescale in wall units
lambda_z = lambda_z / delta_nu

# Print the integral length scale value
print("Integral length scale in viscous units, lambda_z^+ = ", lambda_z)
print()

#!--- Writing to file the viscous time unit and the Kolmogorov time scale ---!
           
# Create the file and write  
with open('data_post/time_scales.txt', 'w') as f:
    f.write(f"{'t_nu':<{pp.c_w}}, "        +
            f"{'min tau_eta':<{pp.c_w}}\n" )  

    f.write(f"{t_nu:{pp.fs}}, "            +
            f"{tau_eta:{pp.fs}}\n"         )      

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

    if pp.iswitch == 0:
    
        # Lee and Moser (2015)
        k = 0.384
        B = 4.27
    
    elif pp.iswitch == 1:
        
        # Cimarelli ('Turbulence' lecture notes)
        k = 0.37
        B = 5.2

# Von Karman law
y_plus_k = np.linspace(5, 180, 175)
u_plus_k = (1.0 / k) * np.log(y_plus_k) + B

#!--------------------------------------------------------------------------------------!
    
#!--- Plot section, mean velocity profile, with selection dipending on the flow case ---!

# Mean velocity profile
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.1
yliminf = 0.0
ylimsup = 25.0

# Mean velocity profile 
ax.scatter(y_plus[:ny], mean_u[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
# TTBL
if itype == 13:
    
    xlimsup = 520.0
                    
# Channel    
elif itype == 3:

    xlimsup = 300.0
    
    # Lee & Moser (2015)
    ax.plot(y_plus_lm, mean_u_lm, color='C1', linestyle='-', linewidth=pp.lw)
    
    # If wall oscillations are present
    if iswitch_wo == 1:
    
        # Touber & Leschziner (2012)
        #ax.plot(y_plus_touber, mean_u_touber, color='C2', linestyle='-.', linewidth=pp.lw)
        
        # Yao et al. (2019)
        ax.scatter(y_plus_umean_yao, mean_u_yao, marker='^', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='k')
    
# Viscous sublayer and log law
ax.plot(y_plus_vsl, u_plus_vsl, color=pp.grey, linestyle='--', linewidth=pp.lw)
ax.plot(y_plus_k, u_plus_k, color=pp.grey, linestyle='--', linewidth=pp.lw)
    
# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$U^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

# Saving the figure
if itype == 13:
    plt.savefig(f'plots/umean-{snap_numb}_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)
elif itype == 3:
    plt.savefig(f'plots/umean_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)
    
# Show the figure
plt.show()

#!--------------------------------------------------------------------------------------!

#!--- Mean spanwise velocity profile ---!

# Mean spanwise velocity profile
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.0
yliminf = min(mean_w)*1.2    
ylimsup = max(mean_w)*1.2

# Spanwise mean velocity profile
ax.scatter(y[:ny], mean_w[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
# TTBL
if itype == 13:
    
    xlimsup = Ly
    
    # Axes labels
    ax.set_xlabel(r'$y/D$',   fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$W/U_w$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Channel
elif itype == 3:
    
    xlimsup = 1.0
    
    # Axes labels
    ax.set_xlabel(r'$y/h$',   fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$W/U_p$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Saving the figure
if itype == 13:
    plt.savefig(f'plots/wmean-{snap_numb}_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)
elif itype == 3:
    plt.savefig(f'plots/wmean_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)
    
# Show the figure
plt.show()

#!--------------------------------------------------------------------------------------!

# <u'u'>
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.1
yliminf = 0.0
ylimsup = 8.0

# <u'u'>
ax.scatter(y_plus[:ny], var_u[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

# TTBL
if itype == 13:
    
    xlimsup = 520.0
           
# Channel    
elif itype == 3:

    xlimsup = 300.0
    
    # Lee & Moser (2015)
    ax.plot(y_plus_lm, var_u_lm, color='C1', linestyle='-', linewidth=pp.lw)
    
    # If wall oscillations are present
    if iswitch_wo == 1:
    
        # Yao et al. (2019)
        ax.scatter(y_plus_uvar_yao, var_u_yao, marker='^', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='k')
    
# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$\langle u^{\prime 2} \rangle^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

# Saving the figure
if itype == 13:
    plt.savefig(f'plots/uvar-{snap_numb}_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)
elif itype == 3:
    plt.savefig(f'plots/uvar_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)

# Show the figure
plt.show()

#!--------------------------------------------------------------------------------------!

# <v'v'>
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.1
yliminf = 0.0
ylimsup = 0.8

# <v'v'>
ax.scatter(y_plus[:ny], var_v[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

# TTBL
if itype == 13:

    xlimsup = 520.0
            
# Channel    
elif itype == 3:

    xlimsup = 300.0
    
    # Lee & Moser (2015)
    ax.plot(y_plus_lm, var_v_lm, color='C1', linestyle='-', linewidth=pp.lw)
    
    # If wall oscillations are present
    if iswitch_wo == 1:
    
        # Yao et al. (2019)
        ax.scatter(y_plus_vvar_yao, var_v_yao, marker='^', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='k')
    
# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$\langle v^{\prime 2} \rangle^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

# Saving the figure
if itype == 13:
    plt.savefig(f'plots/vvar-{snap_numb}_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)
elif itype == 3:
    plt.savefig(f'plots/vvar_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)

# Show the figure
plt.show()

#!--------------------------------------------------------------------------------------!

# <u'v'>
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.1
yliminf = 0.0
ylimsup = 0.8

# <u'v'>
ax.scatter(y_plus[:ny], mean_uv[:ny], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')
    
# TTBL
if itype == 13:

    xlimsup = 520.0
        
    # y-axis label
    ax.set_ylabel(r'$\langle u^{\prime} v^{\prime}\rangle^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
        
# Channel    
elif itype == 3:

    xlimsup = 300.0
        
    # Lee & Moser (2015)
    ax.plot(y_plus_lm, mean_uv_lm, color='C1', linestyle='-', linewidth=pp.lw)
       
    # y-axis label
    ax.set_ylabel(r'$-\langle u^{\prime} v^{\prime}\rangle^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    # If wall oscillations are present
    if iswitch_wo == 1:
    
        # Yao et al. (2019)
        ax.scatter(y_plus_uvmean_yao, mean_uv_yao, marker='^', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='k')
 
# Axes labels
ax.set_xlabel(r'$y^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)

# Saving the figure
if itype == 13:
    plt.savefig(f'plots/uvmean-{snap_numb}_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)
elif itype == 3:
    plt.savefig(f'plots/uvmean_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)

# Show the figure
plt.show()

#!--------------------------------------------------------------------------------------!

# Correlation coefficients in spanwise direction (z)

#!--------------------------------------------------------------------------------------!

# Cuuz
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.0
xlimsup = Lz_plus / 2.0
yliminf = min(Ruuz[c,:])*1.2
ylimsup = max(Ruuz[c,:])*1.2

# Streamwise velocity auto-correlations, Cuuz
ax.scatter(rz, Ruuz[c,:nz], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

# Kim et al. (1987) data
ax.plot(rz_plus_cuuz_kim, cuuz_kim, color='C1', linestyle='-', linewidth=pp.lw)

# Plot horizontal line at Cuu = 0
ax.hlines(y=0.0, xmin=xliminf, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed')
                    
# Axes labels
ax.set_xlabel(r'$r_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$C_{uu}(r_z^+)$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Saving the figure
if itype == 13:
    plt.savefig(f'plots/Cuuz-{snap_numb}_{add_string}_y+={y_plus_in}.pdf', format='pdf', bbox_inches='tight', dpi=600)
elif itype == 3:
    plt.savefig(f'plots/Cuuz_{add_string}_y+={y_plus_in}.pdf', format='pdf', bbox_inches='tight', dpi=600)

# Show the figure
plt.show()

#!--------------------------------------------------------------------------------------!

# Cvvz
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.0
xlimsup = Lz_plus / 2.0
yliminf = min(Rvvz[c,:])*1.2
ylimsup = max(Rvvz[c,:])*1.2

# Vertical velocity auto-correlations, Cvvz
ax.scatter(rz, Rvvz[c,:nz], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

# Kim et al. (1987) data
ax.plot(rz_plus_cvvz_kim, cvvz_kim, color='C1', linestyle='-', linewidth=pp.lw)

# Plot horizontal line at Cvv = 0
ax.hlines(y=0.0, xmin=xliminf, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed')
                  
# Axes labels
ax.set_xlabel(r'$r_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$C_{vv}(r_z^+)$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Saving the figure
if itype == 13:
    plt.savefig(f'plots/Cvvz-{snap_numb}_{add_string}_y+={y_plus_in}.pdf', format='pdf', bbox_inches='tight', dpi=600)
elif itype == 3:
    plt.savefig(f'plots/Cvvz_{add_string}_y+={y_plus_in}.pdf', format='pdf', bbox_inches='tight', dpi=600)

# Show the figure
plt.show()

#!--------------------------------------------------------------------------------------!

# Cwwz
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.0
xlimsup = Lz_plus / 2.0
yliminf = min(Rwwz[c,:])*1.2
ylimsup = max(Rwwz[c,:])*1.2

# Spanwise velocity auto-correlations, Cuuz
ax.scatter(rz, Rwwz[c,:nz], marker='o', linewidth=pp.lw, s=pp.markersize, facecolors='none', edgecolors='C0')

# Kim et al. (1987) data
ax.plot(rz_plus_cwwz_kim, cwwz_kim, color='C1', linestyle='-', linewidth=pp.lw)

# Plot horizontal line at Cww = 0
ax.hlines(y=0.0, xmin=xliminf, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed')
                   
# Axes labels
ax.set_xlabel(r'$r_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$C_{ww}(r_z^+)$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Saving the figure
if itype == 13:
    plt.savefig(f'plots/Cwwz-{snap_numb}_{add_string}_y+={y_plus_in}.pdf', format='pdf', bbox_inches='tight', dpi=600)
elif itype == 3:
    plt.savefig(f'plots/Cwwz_{add_string}_y+={y_plus_in}.pdf', format='pdf', bbox_inches='tight', dpi=600)

# Show the figure
plt.show()

#!--------------------------------------------------------------------------------------!


