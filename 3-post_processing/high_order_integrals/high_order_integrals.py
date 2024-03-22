#!---------------------------------------------------------!
#!    With this script, we perform 6th order accurate      !
#!       calculations of BL thickness parameters           !
#! (delta_99, displacement thickness, momentum thickness), !
#!     related Reynolds numbers and shear velocity.        !
#!---------------------------------------------------------!

# Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

# Parameters
uwall = np.float64(1.0)
re = np.float64(500.0)  
nu = 1/re
ii = 0

# Inputs
t1  = 20  # initial time unit
tn  = 40  # final time unit
icr = 20  # increment between different savings in time units

# Number of snapshots
ns = (tn - t1)//icr + 1 

# Work arrays
delta_99 = np.zeros(ns)
disp_t   = np.zeros(ns)
mom_t    = np.zeros(ns)
sh_vel   = np.zeros(ns)
t_unit   = np.zeros(ns)

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
for i in range(t1, tn + icr, icr):
         
    # Reading of mean streamwise velocity
    file_path = f"mean_stats{i}.0.txt"
    
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
    
    # Reading of the mean gradient
    file_path = f"vort_stats{i}.0.txt"
    
    data = np.loadtxt(file_path, delimiter=',', skiprows=1, dtype=np.float64)
    mg = data[:, 3]
    
    # Shear velocity
    sh_vel[ii] = np.sqrt(nu*np.absolute(mg[0]))
    
    # Time-unit
    t_unit[ii] = t1 + icr*ii
    
    # Index for BL thickness parameters vectors
    ii = ii + 1


# Related Reynolds numbers
re_tau   = delta_99*sh_vel*re
re_ds    = disp_t*uwall*re
re_theta = mom_t*uwall*re

# Combine the vectors into a 2D array
data = np.column_stack((delta_99, disp_t, mom_t, re_tau, re_ds, re_theta, sh_vel, t_unit))

# Specify the file name (ho: High-order) 
filename = 'param_stats_ho.txt'  

# Define the header
header = "delta_99 O(6), disp_t O(6), mom_t O(6), Re_tau O(6), Re_ds O(6), Re_theta O(6), sh_vel O(6), T_unit\n" 

# Write the header and the data to the file
with open(filename, 'w') as file:
     file.write(header)  
     np.savetxt(file, data, delimiter=',')  

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




