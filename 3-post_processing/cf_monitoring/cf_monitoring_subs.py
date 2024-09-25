
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Subroutine(s) used in 'cf_monitoring.py'.  
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

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


"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this subroutine we perform calculation of 6th order 
!              accurate calculations of integrals of TTBL thickness 
!              parameters using 'umean' files printed at the same time as 
!              'cf_monitoring'.  
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def calculate_thickness_param():

    #!--- Reading of files section and setup of flow parameters ---!

    # Read useful flow parameters from 'input.i3d' and 'post.prm' files
    (itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, ioutput, ioutput_cf, iswitch_wo,  
     add_string, file1, filen, icrfile, nr, post_mean, post_vort, post_diss, post_corz, post_tke_eq
    ) = read_input_files('input.i3d','post.prm')

    # Reading of yp coordinates
    yp = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)

    # First element of yp vector (y = 0)
    y0 = yp[0]

    # Last  element of yp vector (y = Ly, height of the domain)    
    yn = yp[-1]

    #!--- Parameters ---!
    uwall, nu, twd = set_flow_parameters(itype, re)
                          
    #!--------------------------------------------------------------------------------------!
                    
    # Number of savings due to 'print_cf' subroutine of Incompact3d solver 'modified'
    nsavings = ilast // ioutput_cf + 1

    # Initialize the mean streamwise velocity profile array (function of y and specific flow realization)
    umean = np.zeros((ny, nr))

    # Initialize the array to sum the mean streamwise velocity profile (function of y and of time)    
    umean_realiz = np.zeros((ny, nsavings))

    # Initialize arrays for TTBL thickness parameters
    disp_t = np.zeros(nsavings)   # displacement thickness, delta*
    mom_t  = np.zeros(nsavings)   # momentum     thickness, theta

    """
    Do loop from the first saving of umean (IC) to the last one, with increment ioutput_cf
    that is read from the input file 'input.i3d'.
    """

    for j in range(0, nsavings, 1):
      
        # Calculate ts to open 'umean-ts' file (ts_iter: time-step of the iterations)
        
        # ts of the initial condition is ts = 1
        if j == 0:

            ts_iter = 1

        # All the other ts are multiples of 'output_cf'
        else:

            ts_iter = j*ioutput_cf
        
        # Do loop over different realizations
        for i in range(1, nr+1, 1):
               
            # Read of 'umean' data from 'data/umean' folder
            umean[:,i-1] = np.loadtxt(f'data_r{i:01d}/umean/umean-ts{ts_iter:07d}.txt', skiprows=1, delimiter=None, dtype=np.float64)
            
            # Summing into a sum array for different realizations
            umean_realiz[:,j] = umean_realiz[:,j] + umean[:,i-1]
    
        #!--- Calculation of thickness parameters ---!

        # Calculate the displacement thickness delta*
        int1 = umean_realiz[:,j]/uwall  # 'integrand 1' 

        # Interpolation at the 6th order of accuracy with a spline of 5th order
        spl = InterpolatedUnivariateSpline(yp, int1, k=5)
        disp_t[j] = spl.integral(y0, yn)

        # Calculate the momentum thickness theta
        int2 = int1 - int1**2  # 'integrand 2' 

        # Interpolation at the 6th order of accuracy with a spline of 5th order
        spl = InterpolatedUnivariateSpline(yp, int2, k=5)
        mom_t[j] = spl.integral(y0, yn)
        
    # Save to file
    

    # Return to main program
    return (disp_t, mom_t)



