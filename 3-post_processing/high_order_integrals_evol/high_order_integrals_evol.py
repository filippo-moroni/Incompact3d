
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this script we perform calculation of 6th order accurate
!              calculation of integrals of TTBL thickness parameters using
!              'umean' files printed at the same time as 'cf_monitoring'.  
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '../../4-common', 'python_common'))
sys.path.append(config_path)

# Import the plotting_params module
import plot_params as pp

# Import functions to setting up, save and show plots a
from plot_subs import set_plot_settings, save_and_show_plot

# Import function to read 'input.i3d' and 'post.prm' files and reference data
from read_files import read_input_files, read_ref_data_temp_evol

# Import function to setup flow parameters 
from set_flow_parameters import set_flow_parameters

#!--------------------------------------------------------------------------------------!

# Print to screen what the program does

print("!--- 'high_order_integrals_evol.py' ---!")
print()
print(" Calculation of:")
print()
print(" TTBL: ")
print(" - displacement thickness, delta*;")
print(" - momentum thickness, theta.")
print()

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
(itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, ioutput, ioutput_cf, iswitch_wo,  
 add_string, file1, filen, icrfile, nr, post_mean, post_vort, post_diss, post_corz, post_tke_eq
) = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!


#!--- Parameters ---!
uwall, nu, twd = set_flow_parameters(itype, re)
                         
#!--- Reading of files section ---!

# Reading of grid points
y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)
    
#!--------------------------------------------------------------------------------------!

# to be completed

# First time-step is skipped at the moment for the reading of umean
ts1 = ioutput_cf
                
# Last time-step
tsn = ilast
    
# Number of savings due to 'print_cf' subroutine of Incompact3d solver 'modified'
nsavings = (tsn - ts1) // ioutput_cf + 1

# Initialize the mean streamwise velocity profile array (function of y and specific flow realization)
umean = np.zeros(ny, nr)

# Initialize the array to sum the mean streamwise velocity profile (function of y and of time)    
umean_realiz = np.zeros(ny, nsavings)
           
"""
Do loop from the first saving of umean (excluding the IC) to the last one, with increment ioutput_cf
that is read from the input file 'input.i3d'.
"""
for j in range(ts1, tsn, ioutput_cf):
        
    # Do loop over different realizations
    for i in range(1, nr + 1, 1):
               
        # Read of 'umean' data from 'data/umean' folder
        umean[:,i] = np.loadtxt(f'data_r{i:01d}/umean/umean-ts{j:07d}.txt', skiprows=1, delimiter=None, dtype=np.float64)
            
        # Summing into a sum array for different realizations
        umean_realiz[:,j] = umean[:,j] + umean[:,i]
    
    # Here, we can calculate the thickness parameters
    
    # ...
    

# Save to file

# Plot (optionally) (we can call this function before cf_monitoring so we can read data later and plot with different available quantities)

print()
print(">>> End.")
print()



