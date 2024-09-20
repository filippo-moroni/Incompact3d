
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this script, we perform 6th order accurate calculations 
!              of thickness parameters (delta_99, displacement thickness, 
!              momentum thickness), related Reynolds numbers, streamwise 
!              shear velocity, streamwise friction coefficient and analogy 
!              factor for a TTBL.
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

print("!--- 'high_order_integrals.py' ---!")
print()
print(" Calculation of:")
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
#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
(itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, ioutput, iswitch_wo,  
 add_string, file1, filen, icrfile, nr, post_mean, post_vort, post_diss, post_corz, post_tke_eq
) = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!

# Create the folder to store the results
os.makedirs('data_post/integral_statistics', mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

#!--- Parameters ---!
uwall, nu, twd = set_flow_parameters(itype, re)

# Local variables
ii    = 0                             # Index for BL thickness parameters vectors 
ns    = (filen - file1)//icrfile + 1  # Number of snapshots

# Work arrays
delta_99  = np.zeros(ns)
disp_t    = np.zeros(ns)
mom_t     = np.zeros(ns)
sh_velx   = np.zeros(ns)
cfx       = np.zeros(ns)
a_fact    = np.zeros(ns)
time_unit = np.zeros(ns)

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
    
#!---------------------------------------------------------!
# Calculations start here, we are employing a Python 
# spline function that passes through all provided points.
#!---------------------------------------------------------!

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
    
    # Index for BL thickness parameters vectors
    ii = ii + 1


# Related Reynolds numbers
re_tau   = delta_99*sh_velx*re
re_ds    = disp_t*uwall*re
re_theta = mom_t*uwall*re

print(">>> Writing to .txt file.")
print()

# Create the file and write  
with open('data_post/integral_statistics/integral_statistics.txt', 'w') as f:
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
            
# Print that calculations have been completed
print(">>> Done!")
print()
print(">>> Results saved in: data_post/integral_statistics/integral_statistics.txt.")
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




