#!---------------------------------------------------------!
#!    With this script, we perform 6th order accurate      !
#!          calculations of thickness parameters           !
#! (delta_99, displacement thickness, momentum thickness), !
#! related Reynolds numbers, streamwise shear velocity and !
#!      streamwise friction coefficient for a TTBL.        !
#!---------------------------------------------------------!

# Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
import os 

# Settings
np.seterr(divide='ignore', invalid='ignore')

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Computer Modern Sans serif",
})

plt.rcParams.update({'figure.autolayout': True})

# Parameters
uwall = np.float64(1.0)
re = np.float64(500.0)  
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
time_unit = np.zeros(ns)

# Reading of yp coordinates
file_path = 'yp.dat'
data = np.loadtxt(file_path, delimiter=None, dtype=np.float64)
yp = data[:]

y0 = data[0]   # First element of yp vector (y = 0)
yn = data[-1]  # Last  element of yp vector (y = Ly, height of the domain)

# Reading of time units (to be fixed: cf is saved more frequently than snapshots!)
# open the file in restart_info folder is a possible solution if checkpoints and snapshots are saved at the same time
#file_path = f"monitoring/cf_history.txt"
      
#data = np.loadtxt(file_path, delimiter=',', skiprows=1, dtype=np.float64)
#time_unit = data[:, 7]

#!---------------------------------------------------------!
# Calculations start here, we are employing a Python 
# spline function that passes through all provided points.
#!---------------------------------------------------------!

# Do loop over different time units
for i in range(file1, filen + icrfile, icrfile):
    
    # Use the first realization folder (named differently if more than 1 realizations are present),
    # to read time units from snapshots
    if nr != 1:
        
        file_path = f"data_r{nr:01d}/snapshot-{i:04d}.xdmf"

    elif nr == 1:
    
        file_path = f"data/snapshot-{i:04d}.xdmf"
    
    # Reading of .xdmf snapshot header
    with open(file_path, 'r') as file:
    
        # Read all lines into a list
        lines = file.readlines()
    
        # Extract time unit line
        input_str    = lines[53]
        
        # Extract the part between the double quotes
        number_str = input_str.split('"')[1].strip()

        # Convert the string to a float
        number = float(number_str)

        print(f"Extracted number: {number}")
    
        # Removing characters around the time unit value          
        #tu    = tu.split('" /')[0]
        #tu    = tu.split('="')[-1].strip()
        #tu    = np.float64(tu)
        
        # Save the extracted time unit value inside its array
        time_unit[ii] = tu 
            
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
    
    # Reading of the mean streamwise gradient
    file_path = f"data_post/vort_stats-{i:04d}.txt"
    
    data = np.loadtxt(file_path, delimiter=',', skiprows=1, dtype=np.float64)
    mg = data[:, 4]
    
    # Shear velocity
    sh_vel[ii] = np.sqrt(nu * (np.absolute(mg[0])))
    
    # Friction coefficient in x-dir.
    cf[ii] = (2.0 * ((sh_vel[ii] / uwall)**2))
            
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




