
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
from scipy.interpolate import InterpolatedUnivariateSpline

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '../../4-common', 'python_common'))
sys.path.append(config_path)

# Import function to read 'input.i3d' and 'post.prm' files
from read_files import read_input_files

# Import function to setup flow parameters 
from set_flow_parameters import set_flow_parameters

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this subroutine we perform 6th order accurate 
!              calculations of integrals of TTBL thickness parameters 
!              (delta*, theta) using 'umean' files printed at the same time 
!              as 'cf_monitoring'.  
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

    # Initialize instantaneous mean statistics array, not averaged with different flow realizations
    # we have 9 different statistics (mean, var, Reynolds stress) (function of y)
    mean_stats = np.zeros((ny, 9))
    
    # Initialize mean statistics array, function of y and time
    mean_stats_realiz = np.zeros((ny, 9, nsavings))

    # Initialize arrays for TTBL thickness parameters
    disp_t = np.zeros(nsavings)   # displacement thickness, delta*
    mom_t  = np.zeros(nsavings)   # momentum     thickness, theta

    """
    Do loop from the first saving of umean (IC) to the last one, with increment ioutput_cf.
    """
    for j in range(0, nsavings, 1):
      
        #!--- Calculate ts to open 'umean-ts' file (ts_iter: time-step of the iterations) ---!
        
        # ts of the initial condition is ts = 1
        if j == 0:

            ts_iter = 1

        # All the other ts are multiples of 'output_cf'
        else:

            ts_iter = j*ioutput_cf
        
        # Do loop over different realizations, from 1 to nr
        for i in range(1, nr+1, 1):
               
            # Read of mean statistics calculated runtime data from 'data/mean_stats_runtime' folder
            mean_stats = np.loadtxt(f'data_r{i:01d}/mean_stats_runtime/mean_stats_runtime-ts{ts_iter:07d}.txt', skiprows=10, delimiter=None, dtype=np.float64)
                        
            # Summing mean statistics array with different realizations into the overall array for time-evolution
            mean_stats_realiz[:,:,j] = mean_stats_realiz[:,:,j] + mean_stats[:,:] / nr
            
        #!--- Finalize 2nd order statistics ---!
        
        # Variances
        mean_stats_realiz[:,3,:] = mean_stats_realiz[:,3,:] - mean_stats_realiz[:,0,:]**2  # streamwise  velocity variance
        mean_stats_realiz[:,4,:] = mean_stats_realiz[:,4,:] - mean_stats_realiz[:,1,:]**2  # wall-normal velocity variance
        mean_stats_realiz[:,5,:] = mean_stats_realiz[:,5,:] - mean_stats_realiz[:,2,:]**2  # spanwise    velocity variance    
        
        # Reynolds stress
        mean_stats_realiz[:,6,:] = mean_stats_realiz[:,6,:] - mean_stats_realiz[:,0,:]*mean_stats_realiz[:,1,:]  # Reynolds stress <u'v'>
        mean_stats_realiz[:,7,:] = mean_stats_realiz[:,7,:] - mean_stats_realiz[:,0,:]*mean_stats_realiz[:,2,:]  # Reynolds stress <u'w'>
        mean_stats_realiz[:,8,:] = mean_stats_realiz[:,8,:] - mean_stats_realiz[:,1,:]*mean_stats_realiz[:,2,:]  # Reynolds stress <v'w'>

        #!--- Calculation of thickness parameters ---!

        # Calculate the displacement thickness delta*
        # First column of the 'mean_stats_realiz' array is mean streamwise velocity profile
        int1 = mean_stats_realiz[:,0,j]/uwall  # 'integrand 1' 

        # Interpolation at the 6th order of accuracy with a spline of 5th order, 
        # that passes through all data points
        spl = InterpolatedUnivariateSpline(yp, int1, k=5)
        disp_t[j] = spl.integral(y0, yn)

        # Calculate the momentum thickness theta
        int2 = int1 - int1**2  # 'integrand 2' 

        # Interpolation at the 6th order of accuracy with a spline of 5th order,
        # that passes through all data points
        spl = InterpolatedUnivariateSpline(yp, int2, k=5)
        mom_t[j] = spl.integral(y0, yn)

    # Return to main program
    return (disp_t, mom_t)



