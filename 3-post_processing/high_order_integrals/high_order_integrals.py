#!---------------------------------------------------------!
#!    With this script, we perform 6th order accurate      !
#!          calculations of thickness parameters           !
#! (delta_99, displacement thickness, momentum thickness), !
#! related Reynolds numbers, streamwise shear velocity,    !
#!  streamwise friction coefficient and analogy factor     ! 
#!                       for a TTBL.                       !
#!---------------------------------------------------------!

# We are assuming unitary molecular Prandtl number for the calculation
# of the analogy factor of the Reynolds analogy. 

# Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
import os
import xml.etree.ElementTree as ET 

# Settings
np.seterr(divide='ignore', invalid='ignore')

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Computer Modern Sans serif",
})

plt.rcParams.update({'figure.autolayout': True})

# Parameters
uwall = np.float64(1.0)    # Wall velocity, Uw
re    = np.float64(500.0)  # Trip Reynolds number, Re_D
nu = 1.0/re
ii = 0

# Reading of the post-processing input file for indexes calculation of snapshots
file_path = f"post.prm"
data = np.loadtxt(file_path, delimiter=None, dtype=int, comments='#', skiprows=5)

# Inputs
file1   = data[0]  # First snapshot index
filen   = data[1]  # Final snapshot index
icrfile = data[2]  # File increment
nr      = data[3]  # Number of flow realizations

# Number of snapshots
ns = (filen - file1)//icrfile + 1 

# Work arrays
delta_99  = np.zeros(ns)
disp_t    = np.zeros(ns)
mom_t     = np.zeros(ns)
sh_vel    = np.zeros(ns)
cf        = np.zeros(ns)
a_fact    = np.zeros(ns)
time_unit = np.zeros(ns)

# Reading of yp coordinates
file_path = 'yp.dat'
data = np.loadtxt(file_path, delimiter=None, dtype=np.float64)
yp = data[:]

y0 = data[0]   # First element of yp vector (y = 0)
yn = data[-1]  # Last  element of yp vector (y = Ly, height of the domain)

#!---------------------------------------------------------!
# Calculations start here, we are employing a Python 
# spline function that passes through all provided points.
#!---------------------------------------------------------!

# Do loop over different time units
for i in range(file1, filen + icrfile, icrfile):
    
    # Use first realization folder to read time units from snapshots
    if nr != 1:
        # /data_ri folders if nr /= 1
        file_path = f"data_r{nr:01d}/snapshot-{i:04d}.xdmf"
       
    elif nr == 1:
        # /data folder if nr = 1
        file_path = f"data/snapshot-{i:04d}.xdmf"
    
    # Read .xdmf file of snapshots to extract time unit
    tree = ET.parse(file_path)
    root = tree.getroot()
        
    # Find the 'Time' attribute within the 'Grid' element
    time_element = root.find(".//Grid/Time")
    if time_element is not None:
        time_str = time_element.attrib.get('Value', '').strip()
        try:
            time_unit[ii] = float(time_str)
            print(f"Extracted time unit for snapshot-{i:04d}: {time_unit[ii]}")
        except ValueError:
            print(f"Warning: Failed to convert '{time_str}' to float for snapshot-{i:04d}")
    else:
        print(f"Warning: 'Time' element not found in {file_path}")
                
    # Reading of mean streamwise velocity
    file_path = f"data_post/mean_stats-{i:04d}.txt"
       
    data = np.loadtxt(file_path, delimiter=',', skiprows=1, dtype=np.float64)
    umean = data[:, 0]
    
    # Calculate the BL thickness delta99
    j = 0
    while umean[j] > umean[0]*0.01:
        delta_99[ii] = yp[j]
        j = j + 1 

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
        
    # Shear velocity
    sh_vel[ii] = np.sqrt(nu * (np.absolute(mgx[0])))
    
    # Friction coefficient in x-dir.
    cf[ii] = (2.0 * ((sh_vel[ii] / uwall)**2))
            
    # Analogy factor, ratio between mean gradient parallel to the wall of velocity and mean scalar gradient
    a_fact[ii] = mgpar[0] / mgphi[0]
    
    # Index for BL thickness parameters vectors
    ii = ii + 1


# Related Reynolds numbers
re_tau   = delta_99*sh_vel*re
re_ds    = disp_t*uwall*re
re_theta = mom_t*uwall*re

# Column width for writing to .txt file
c_w = 20  

# Format for numbers
fs = f"<{c_w}.3f"

# Format for cf only
fs2 = f"<{c_w}.6f"

# Create the folder to store the results
os.makedirs('integral_statistics', mode=0o777, exist_ok=True)

# Create the file and write  
with open('integral_statistics/integral_statistics.txt', 'w') as f:
    f.write(f"{'delta_99 O(6)':<{c_w}}, " +
            f"{'disp_t O(6)':<{c_w}}, " +
            f"{'mom_t O(6)':<{c_w}}, " +
            f"{'Re_tau O(6)':<{c_w}}, " +
            f"{'Re_ds O(6)':<{c_w}}, " +
            f"{'Re_theta O(6)':<{c_w}}, " +
            f"{'sh_velx O(6)':<{c_w}}, " +
            f"{'cf,x O(6)':<{c_w}}, " +
            f"{'A_fact O(6)':<{c_w}}, " +
            f"{'time_unit':<{c_w}}\n")

    for j in range(0, ii):
        f.write(f"{delta_99[j]:{fs}}, " +
            f"{disp_t[j]:{fs}}, " +
            f"{mom_t[j]:{fs}}, " +
            f"{re_tau[j]:{fs}}, " +
            f"{re_ds[j]:{fs}}, " +
            f"{re_theta[j]:{fs}}, " +
            f"{sh_vel[j]:{fs}}, " +
            f"{cf[j]:{fs2}}, " +
            f"{a_fact[j]:{fs}}, " +
            f"{time_unit[j]:{fs}}\n")
           
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




