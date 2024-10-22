
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this script, we perform 6th order accurate calculations 
!              of thickness parameters (delta_99, displacement thickness, 
!              momentum thickness), related Reynolds numbers, streamwise 
!              shear velocity, streamwise friction coefficient and analogy 
!              factor for a TTBL.
!              We are also calculating non-dimensional grid spacings and 
!              domain dimensions at each snapshots' saving. Finally, we 
!              calculate minimum Kolmogorov time scale and viscous time unit.
! ANNOTATIONS: 1) We are assuming unitary molecular Prandtl number for the 
!                 calculation of the analogy factor of the Reynolds analogy.
!              2) We are calculating these quantities through averages of
!                 snapshots. For a complete temporal evolution, 
!                 use 'cf_monitoring.py'.                         
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

# Import function to calculate TTBL thickness parameters
from ttbl_subs import calculate_ttbl_thick_params

#!--------------------------------------------------------------------------------------!

# Print to screen what the program does

print("!--- 'ttbl_indexes_snap.py' ---!")
print()
print(" File 'data_post/ttbl_indexes_snap/thickness_params_evolution.txt':")
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
print(" File 'data_post/ttbl_indexes_snap/nd_mesh_evolution.txt':")
print("  - non-dimensional mesh spacings;")
print("  - non-dimensional domain dimensions.")
print()
print()
print(" File 'data_post/ttbl_indexes_snap/time_scales_evolution.txt':")
print("  - minimum Kolmogorov time scale;")
print("  - viscous time unit.")
print()
print(" The calculated quantities are checked at each snapshot only.")
print(" For a complete temporal evolution, use 'cf_monitoring.py'.")
print()

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
(itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, itimescheme, ioutput, ioutput_cf, iswitch_wo,  
 add_string, file1, filen, icrfile, nr, post_mean, post_grad, post_corz, post_tke_eq
) = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!

# Create the folder to store the results
os.makedirs('time_evolution',  mode=0o777, exist_ok=True)
os.makedirs('num_resolutions', mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

#!--- Parameters and mesh ---!
(nu, y, uwall, twd, phiwall) = set_flow_parameters(re)

# Mesh spacings
delta_x = Lx / nx
delta_z = Lz / nz

# Local variables
ii    = 0                             # Index for BL thickness parameters vectors 
ns    = (filen - file1)//icrfile + 1  # Number of snapshots

# TTBL thickness parameters and shear quantities
delta_99   = np.zeros(ns)
delta_99_j = np.zeros(ns)
disp_t     = np.zeros(ns)
mom_t      = np.zeros(ns)
sh_velx    = np.zeros(ns)
sh_vel_tot = np.zeros(ns)
cfx        = np.zeros(ns)
a_fact     = np.zeros(ns)
time_unit  = np.zeros(ns)

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
!------------------------------------------------------------------!
 Calculations start here, for the calculation of TTBL thickness
 parameters, we are employing a SciPy spline function that passes 
 through all provided points.
!------------------------------------------------------------------!
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
    mean_u = data[:, 0]
    
    # Call of external subroutine for the calculation of TTBL thickness parameters
    (delta_99[ii], dummy, disp_t[ii], mom_t[ii]) = calculate_ttbl_thick_params(mean_u,y,uwall)
            
    # Reading of mean gradients: streamwise, spanwise and scalar 
    file_path = f"data_post/vort_stats-{i:04d}.txt"
    
    data = np.loadtxt(file_path, delimiter=',', skiprows=1, dtype=np.float64)
    
    mgx   = data[:, 3]  # mean streamwise gradient
    mgz   = data[:, 4]  # mean spanwise   gradient
    mgphi = data[:, 5]  # mean scalar     gradient

    # Mean gradient parallel to the wall
    mg_parw = np.sqrt(mgx[0]**2 + mgz[0]**2)
        
    # Streamwise shear velocity
    sh_velx[ii] = np.sqrt(nu * (np.absolute(mgx[0])))
    
    # Streamwise friction coefficient
    cfx[ii] = (2.0 * ((sh_velx[ii] / uwall)**2))
            
    # Analogy factor, ratio between mean gradient parallel to the wall of velocity and mean scalar gradient
    a_fact[ii] = np.abs(mgphi[0] / mg_parw)
    
    
    """
    !------------------------------------------------------------------!
     Extra section for calculations of grid spacings.
    !------------------------------------------------------------------!    
    """
    
    # Total shear velocity
    sh_vel_tot[ii] = np.sqrt(nu * (np.absolute(mg_parw)))
    
    # Viscous length
    delta_nu = nu / sh_vel_tot[ii]
    
    # Rescaling variables through wall units
    delta_x_plus[ii] = delta_x / delta_nu
    delta_z_plus[ii] = delta_z / delta_nu
    
    y_plus           = y       / delta_nu
    
    delta_yw_plus[ii] = y_plus[1]
    
    # Delta y+ at the BL edge
    delta_yd_plus[ii] = y_plus[bl_thick_j] - y_plus[bl_thick_j-1] 
 
    Lx_plus[ii] = Lx / delta_nu
    Ly_plus[ii] = Ly / delta_nu 
    Lz_plus[ii] = Lz / delta_nu
    
    
    """
    !------------------------------------------------------------------!
     Extra section to calculate time scales: 
     minimum Kolmogorov time scale and viscous time unit.
    !------------------------------------------------------------------!
    """
    
    # Reading of total dissipation
    file_path = f"data_post/diss_stats-{i:04d}.txt"
    
    eps = np.loadtxt(file_path, delimiter=',', skiprows=1, dtype=np.float64)
        
    # Find the maximum of mean total dissipation
    eps_max = max(eps)

    # Minimum Kolmogorov time scale
    tau_eta[ii] = np.sqrt(nu/eps_max)
    
    # Viscous time unit
    t_nu[ii] = nu / (sh_vel_tot[ii] ** 2)
    
    #!-----!
    
    # Index to advance in time along different snapshots
    ii = ii + 1
    
                   
# Related Reynolds numbers
re_tau   = delta_99*sh_velx*re
re_ds    = disp_t*uwall*re
re_theta = mom_t*uwall*re

#!--- Create the files and write ---!

print(">>> Writing to .txt file.")
print()

# Integral statistics and flow indexes at snapshots' time of saving
with open('time_evolution/time_evolution_snaps.txt', 'w') as f:
    f.write('Time evolution of Temporal Turbulent Boundary Layer (TTBL) thickness parameters, \n')    
    f.write('shear quantities and Reynolds analogy factor at each snapshot.\n')        
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
    f.write(f"{'delta_99':>{pp.c_w}}, "     +
            f"{'disp_t':>{pp.c_w}}, "       +
            f"{'mom_t':>{pp.c_w}}, "        +
            f"{'Re_tau':>{pp.c_w}}, "       +    
            f"{'Re_delta*':>{pp.c_w}}, "    +
            f"{'Re_theta:>{pp.c_w}}, "      +
            f"{'sh_velx':>{pp.c_w}}, "      +
            f"{'cfx':>{pp.c_w}}, "          +
            f"{'A_fact:>{pp.c_w}}, "        +
            f"{'time_unit':>{pp.c_w}}\n"    )

    for j in range(0, ii):
        f.write(f"{delta_99[j]:{pp.fs}}, "  +
                f"{disp_t[j]:{pp.fs}}, "    +
                f"{mom_t[j]:{pp.fs}}, "     +
                f"{re_tau[j]:{pp.fs}}, "    +
                f"{re_ds[j]:{pp.fs}}, "     +
                f"{re_theta[j]:{pp.fs}}, "  +
                f"{sh_velx[j]:{pp.fs6}}, "  +
                f"{cfx[j]:{pp.fs8}}, "      +
                f"{a_fact[j]:{pp.fs}}, "    +
                f"{time_unit[j]:{pp.fs}}\n" )
                
# Non-dimensional grid spacings and domain dimensions (nd: non-dimensional)
with open('num_resolutions/nd_mesh_evolution.txt', 'w') as f:
    f.write('Time evolution of Temporal Turbulent Boundary Layer (TTBL) non-dimensional\n')    
    f.write('grid spacings and domain dimensions at each snapshot.\n')
    f.write('Adimensionalization in viscous units (^+) with the (total) shear velocity.\n')        
    f.write('\n')
    f.write(f'Flowcase: {add_string}.\n')
    f.write('\n') 
    f.write('Abbreviations:\n')
    f.write(' - x:        streamwise direction;\n')
    f.write(' - y:        wall-normal direction;\n')
    f.write(' - z:        spanwise direction;\n')
    f.write(' - delta:    mesh spacing;\n')
    f.write(' - L:        domain dimension;\n')
    f.write(' - d:        boundary layer interface (d: small letter greek delta);\n')
    f.write(' - Re_tau:   (streamwise) friction Reynolds number.\n')
    f.write('\n')
    f.write(f"{'delta_x^+':>{pp.c_w}}, "     +
            f"{'delta_yw^+':>{pp.c_w}}, "    +
            f"{'delta_yd^+':>{pp.c_w}}, "    +
            f"{'delta_z^+':>{pp.c_w}}, "     +
            f"{'Lx^+':>{pp.c_w}}, "          +
            f"{'Ly^+':>{pp.c_w}}, "          +
            f"{'Lz^+':>{pp.c_w}}, "          +
            f"{'Re_tau O(6)':>{pp.c_w}}\n"   )  

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
with open('num_resolutions/time_scales_evolution.txt', 'w') as f:
    f.write('Time evolution of Temporal Turbulent Boundary Layer (TTBL) minimum\n')    
    f.write('Kolmogorov time scale and viscous time unit at each snapshot.\n')        
    f.write('\n')
    f.write(f'Flowcase: {add_string}.\n')
    f.write('\n') 
    f.write('Abbreviations:\n')
    f.write(' - tau_eta:  (minimum) Kolmogorov time scale;\n')
    f.write(' - t_nu:     viscous time unit (based on total shear velocity).\n')
    f.write(' - Re_tau:   (streamwise) friction Reynolds number.\n')
    f.write('\n')
    f.write(f'For reference, time-step dt = {dt}.\n')
    f.write('\n')
    f.write(f"{'tau_eta':>{pp.c_w}}, "       +
            f"{'t_nu':>{pp.c_w}}, "          +
            f"{'Re_tau O(6)':>{pp.c_w}}\n"   )  

    for j in range(0, ii):
        f.write(f"{tau_eta[j]:{pp.fs}}, "  +
                f"{t_nu[j]:{pp.fs}}, "     +
                f"{re_tau[j]:{pp.fs}}\n"   )

# Print that calculations have been completed
print(">>> Done!")
print()
print(">>> Results saved in /time_evolution and /num_resolutions.")
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




