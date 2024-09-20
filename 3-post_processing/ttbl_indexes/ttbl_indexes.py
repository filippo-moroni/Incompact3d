
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this script, we perform 6th order accurate calculations 
!              of thickness parameters (delta_99, displacement thickness, 
!              momentum thickness), related Reynolds numbers, streamwise 
!              shear velocity, streamwise friction coefficient and analogy 
!              factor for a TTBL.
!              We are also calculating non-dimensional grid spacings and 
!              domain dimensions at each snapshots' saving.
! ANNOTATIONS: We are assuming unitary molecular Prandtl number for the 
!              calculation of the analogy factor of the Reynolds analogy.                        
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET 
from scipy.interpolate import InterpolatedUnivariateSpline

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '../../4-common', 'python_common'))
sys.path.append(config_path)

# Import the plotting_params module
import plot_params as pp

# Import function to read 'input.i3d' and 'post.prm' files
from read_files import read_input_files

# Import function to setup flow parameters 
from set_flow_parameters import set_flow_parameters

# Import function to calculate boundary layer thickness delta_99 for a TTBL
from ttbl_subs import calculate_ttbl_delta_99

#!--------------------------------------------------------------------------------------!

# Print to screen what the program does

print("!--- 'ttbl_indexes.py' ---!")
print()
print(" File 'data_post/ttbl_indexes/thickness_params_evolution.txt':")
print("  - delta_99;")
print("  - displacement thickness, delta*;")
print("  - momentum thickness, theta;")
print("  - related Re numbers;")
print("  - (streamwise) shear velocity, sh_velx;")
print("  - (streamwise) friction coefficient, cfx;")
print("  - Reynolds analogy factor, A_fact.")
print()
print(" The mean velocity profile is interpolated with a spline of order 5")
print(" that is constrained to the calculated values.")
print()
print(" This function requires at least one realization folder with")
print(" snapshots' headers in order to read time unit t.")
print()
print()
print(" File 'data_post/ttbl_indexes/nd_mesh_evolution.txt':")
print("  - non-dimensional mesh spacings;")
print("  - non-dimensional domain dimensions.")
print()
print()
print(" File 'data_post/ttbl_indexes/time_scales_evolution.txt':")
print("  - minimum Kolmogorov time scale;")
print("  - viscous time unit.")
print()

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
(itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, ioutput, iswitch_wo,  
 add_string, file1, filen, icrfile, nr, post_mean, post_vort, post_diss, post_corz, post_tke_eq
) = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!

# Create the folder to store the results
os.makedirs('data_post/ttbl_indexes', mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

#!--- Parameters ---!
uwall, nu, twd = set_flow_parameters(itype, re)

# Mesh spacings
delta_x = Lx / nx
delta_z = Lz / nz

# Local variables
ii    = 0                             # Index for BL thickness parameters vectors 
ns    = (filen - file1)//icrfile + 1  # Number of snapshots

# TTBL thickness parameters and shear quantities
delta_99  = np.zeros(ns)
disp_t    = np.zeros(ns)
mom_t     = np.zeros(ns)
sh_velx   = np.zeros(ns)
cfx       = np.zeros(ns)
a_fact    = np.zeros(ns)
time_unit = np.zeros(ns)

# Arrays for non-dimensional grid spacings and domain dimensions
delta_x_plus  = np.zeros(ns)
delta_yw_plus = np.zeros(ns)
delta_yd_plus = np.zeros(ns)
delta_z_plus  = np.zeros(ns)
Lx_plus       = np.zeros(ns)
Ly_plus       = np.zeros(ns)
Lz_plus       = np.zeros(ns)

# Arrays for Kolmogorov time scale and viscous time unit
tau_eta       = np.zeros(ns)
t_nu          = np.zeros(ns)

# Reading of yp coordinates
yp = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)
y0 = yp[0]   # First element of yp vector (y = 0)
yn = yp[-1]  # Last  element of yp vector (y = Ly, height of the domain)

#!--------------------------------------------------------------------------------------!

# Path for generic data
data_path = 'data'

# Check if the path exists and is a directory
if os.path.exists(data_path) and os.path.isdir(data_path):
    
    # Confirm /data as path for reading .xdmf files
    data_path = 'data'
        
else:

    # Asking the user the realization folder to use (if TTBL)
    print()
    realiz = int(input(">>> Specify the realization folder to read time units 't': "))
    print()
    data_path = f'data_r{realiz}'

    
"""
Calculations start here, for the calculation of TTBL thickness
parameters, we are employing a SciPy spline function that passes 
through all provided points.

"""

print(">>> Calculations start now.")
print()

# Do loop over different time units
for i in range(file1, filen + icrfile, icrfile):
    
    # Create dynamic file path for .xdmf files
    file_path = data_path + f"/snapshot-{i:04d}.xdmf"
    
    # Read .xdmf file of snapshots to extract time unit
    tree = ET.parse(file_path)
    root = tree.getroot()
        
    # Find the 'Time' attribute within the 'Grid' element
    time_element = root.find(".//Grid/Time")
    if time_element is not None:
        time_str = time_element.attrib.get('Value', '').strip()
        time_unit[ii] = np.float64(time_str)
                        
    # Reading of mean streamwise velocity
    file_path = f"data_post/mean_stats-{i:04d}.txt"
       
    data = np.loadtxt(file_path, delimiter=',', skiprows=1, dtype=np.float64)
    umean = data[:, 0]
    
    # Calculate BL thickness delta_99 for a TTBL and its related index
    (bl_thick, bl_thick_j) = calculate_ttbl_delta_99(umean, yp) 
    
    # Copy the result in the temporal array
    delta_99[ii] = bl_thick

    # Calculate the displacement thickness delta*
    int1 = umean/uwall  # 'integrand 1' 

    # Interpolation at the 6th order of accuracy with a spline of 5th order
    spl = InterpolatedUnivariateSpline(yp, int1, k=5)
    disp_t[ii] = spl.integral(y0, yn)

    # Calculate the momentum thickness theta
    int2 = int1 - int1**2  # 'integrand 2' 

    # Interpolation at the 6th order of accuracy with a spline of 5th order
    spl = InterpolatedUnivariateSpline(yp, int2, k=5)
    mom_t[ii] = spl.integral(y0, yn)
    
    # Reading of mean parallel gradient, mean streamwise gradient and mean scalar gradient
    file_path = f"data_post/vort_stats-{i:04d}.txt"
    
    data = np.loadtxt(file_path, delimiter=',', skiprows=1, dtype=np.float64)
    
    mgpar = data[:, 3]  # mean gradient parallel to the wall
    mgx   = data[:, 4]  # mean streamwise gradient 
    mgphi = data[:, 6]  # mean scalar gradient
        
    # Streamwise shear velocity
    sh_velx[ii] = np.sqrt(nu * (np.absolute(mgx[0])))
    
    # Streamwise friction coefficient
    cfx[ii] = (2.0 * ((sh_velx[ii] / uwall)**2))
            
    # Analogy factor, ratio between mean gradient parallel to the wall of velocity and mean scalar gradient
    a_fact[ii] = np.abs(mgphi[0] / mgpar[0])
    
    
    """
    Extra section for calculations of grid spacings.
        
    """
    
    # Viscous length
    delta_nu = nu / sh_velx[ii]
    
    # Rescaling variables through wall units
    delta_x_plus[ii] = delta_x / delta_nu
    delta_z_plus[ii] = delta_z / delta_nu
    
    y_plus           = yp      / delta_nu
    
    delta_yw_plus[ii] = y_plus[1]
    
    # Delta y+ at the BL edge
    delta_yd_plus[ii] = y_plus[bl_thick_j] - y_plus[bl_thick_j-1] 
 
    Lx_plus[ii] = Lx / delta_nu
    Ly_plus[ii] = Ly / delta_nu 
    Lz_plus[ii] = Lz / delta_nu
    
    
    """
    Extra section to calculate time scales: 
    minimum Kolmogorov time scale and viscous time unit.
    
    """
    
    # Reading of total dissipation
    file_path = f"data_post/diss_stats-{i:04d}.txt"
    
    eps = np.loadtxt(file_path, delimiter=',', skiprows=1, dtype=np.float64)
        
    # Find the maximum of mean total dissipation
    eps_max = max(eps)

    # Minimum Kolmogorov time scale
    tau_eta[ii] = np.sqrt(nu/eps_max)
    
    # Viscous time unit
    t_nu[ii] = nu / (sh_velx[ii] ** 2)
    
    
    # Index to advance in time along different snapshots
    ii = ii + 1
    
                   
# Related Reynolds numbers
re_tau   = delta_99*sh_velx*re
re_ds    = disp_t*uwall*re
re_theta = mom_t*uwall*re

print(">>> Writing to .txt file.")
print()

#!--- Create the files and write ---!

# Integral statistics and flow indexes
with open('data_post/ttbl_indexes/thickness_params_evolution.txt', 'w') as f:
    f.write('Time evolution of Temporal Turbulent Boundary Layer (TTBL) thickness parameters, \n')    
    f.write('shear quantities and Reynolds analogy factor.\n')        
    f.write('\n')
    f.write(f'Flowcase: {add_string}.\n')
    f.write('\n')    
    f.write('Quantities in this file:\n')
    f.write(' - thickness delta_99:              delta_99;\n')
    f.write(' - displacement thickness:          disp_t;\n')
    f.write(' - momentum thickness:              mom_t;\n')
    f.write(' - related Re numbers;\n')
    f.write(' - streamwise shear velocity:       sh_velx;\n')
    f.write(' - streamwise friction coefficient: cfx;\n')
    f.write(' - Reynolds analogy factor:         A_fact;\n')
    f.write(' - (outer) time unit:               time_unit.\n')
    f.write('\n')
    f.write(f"{'delta_99 O(6)':>{pp.c_w}}, " +
            f"{'disp_t O(6)':>{pp.c_w}}, "   +
            f"{'mom_t O(6)':>{pp.c_w}}, "    +
            f"{'Re_tau O(6)':>{pp.c_w}}, "   +
            f"{'Re_ds O(6)':>{pp.c_w}}, "    +
            f"{'Re_theta O(6)':>{pp.c_w}}, " +
            f"{'sh_velx O(6)':>{pp.c_w}}, "  +
            f"{'cfx O(6)':>{pp.c_w}}, "      +
            f"{'A_fact O(6)':>{pp.c_w}}, "   +
            f"{'time_unit':>{pp.c_w}}\n"     )

    for j in range(0, ii):
        f.write(f"{delta_99[j]:{pp.fs}}, "   +
                f"{disp_t[j]:{pp.fs}}, "     +
                f"{mom_t[j]:{pp.fs}}, "      +
                f"{re_tau[j]:{pp.fs}}, "     +
                f"{re_ds[j]:{pp.fs}}, "      +
                f"{re_theta[j]:{pp.fs}}, "   +
                f"{sh_velx[j]:{pp.fs6}}, "   +
                f"{cfx[j]:{pp.fs8}}, "       +
                f"{a_fact[j]:{pp.fs}}, "     +
                f"{time_unit[j]:{pp.fs}}\n"  )
                
# Non-dimensional grid spacings and domain dimensions (nd: non-dimensional)
with open('data_post/ttbl_indexes/nd_mesh_evolution.txt', 'w') as f:
    f.write('Time evolution of Temporal Turbulent Boundary Layer (TTBL) non-dimensional\n')    
    f.write('grid spacings and domain dimensions. Adimensionalization in viscous units (^+).\n')        
    f.write('\n')
    f.write(f'Flowcase: {add_string}.\n')
    f.write('\n') 
    f.write('Abbreviations:\n')
    f.write(' - x:      streamwise direction;\n')
    f.write(' - y:      wall-normal direction;\n')
    f.write(' - z:      spanwise direction;\n')
    f.write(' - delta:  mesh spacing;\n')
    f.write(' - L:      domain dimension;\n')
    f.write(' - d:      boundary layer interface (d: small letter greek delta);\n')
    f.write(' - Re_tau: friction Reynolds number.\n')
    f.write('\n')
    f.write(f"{'delta_x^+':>{pp.c_w}}, "   +
            f"{'delta_yw^+':>{pp.c_w}}, "  +
            f"{'delta_yd^+':>{pp.c_w}}, "  +
            f"{'delta_z^+':>{pp.c_w}}, "   +
            f"{'Lx^+':>{pp.c_w}}, "        +
            f"{'Ly^+':>{pp.c_w}}, "        +
            f"{'Lz^+':>{pp.c_w}}, "        +
            f"{'Re_tau O(6)':>{pp.c_w}}\n" )  

    for j in range(0, ii):
        f.write(f"{delta_x_plus[j]:{pp.fs}}, "  +
                f"{delta_yw_plus[j]:{pp.fs}}, " +
                f"{delta_yd_plus[j]:{pp.fs}}, " +
                f"{delta_z_plus[j]:{pp.fs}}, "  +
                f"{Lx_plus[j]:{pp.fs}}, "       +
                f"{Ly_plus[j]:{pp.fs}}, "       +
                f"{Lz_plus[j]:{pp.fs}}, "       +
                f"{re_tau[j]:{pp.fs}}\n"        )
                
# Time scales (minimum Kolmogorov time scale and viscous time unit)
with open('data_post/ttbl_indexes/time_scales_evolution.txt', 'w') as f:
    f.write('Time evolution of Temporal Turbulent Boundary Layer (TTBL) minimum\n')    
    f.write('Kolmogorov time scale and viscous time unit.\n')        
    f.write('\n')
    f.write(f'Flowcase: {add_string}.\n')
    f.write('\n') 
    f.write('Abbreviations:\n')
    f.write(' - tau_eta: (minimum) Kolmogorov time scale;\n')
    f.write(' - t_nu:     viscous time unit.\n')
    f.write(' - Re_tau:   friction Reynolds number.\n')
    f.write('\n')
    f.write(f'For reference, time-step dt = {dt}.\n')
    f.write('\n')
    f.write(f"{'tau_eta':>{pp.c_w}}, "     +
            f"{'t_nu':>{pp.c_w}}, "        +
            f"{'Re_tau O(6)':>{pp.c_w}}\n" )  

    for j in range(0, ii):
        f.write(f"{tau_eta[j]:{pp.fs}}, " +
                f"{t_nu[j]:{pp.fs}}, "    +
                f"{re_tau[j]:{pp.fs}}\n"  )

            
# Print that calculations have been completed
print(">>> Done!")
print()
print(">>> Results saved in: data_post/ttbl_indexes.")
print()
           
#!---------------------------------------------------------!
# Check section (everything works)
#!---------------------------------------------------------!

# Very similar values to low-order calculations
#print("Displacement thickness:", disp_t)
#print("Momentum thickness:", mom_t)

# Spline reconstructs well the mean velocity profile
#plt.plot(yp, spl(yp), 'g', lw=3, alpha=0.7, label='spline')
#plt.plot(yp, umean, 'ro', ms=5, label='original data')
#plt.legend()
#plt.show()




