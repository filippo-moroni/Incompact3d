
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

# Import the plotting_params module
import plot_params as pp

# Import function to read 'input.i3d'
from read_files import read_input_files

# Import function to setup flow parameters 
from set_flow_parameters import set_flow_parameters

# Import function to calculate boundary layer thickness delta_99 for a TTBL
from ttbl_subs import calculate_ttbl_delta_99

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this subroutine we perform:
!               - 6th order accurate calculations of integrals of TTBL 
!                 thickness parameters (delta*, theta) using 
!                 'mean_stats_runtime' files printed at the same time as 
!                 'cf_monitoring';
!               - averaging and saving of runtime mean statistics
!                 (different flow realizations) and save related 
!                 shear velocities; 
!               - calculation of TTBL thickness delta_99.
!               
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def calculate_thickness_param(sh_veltot,sh_velx):

    # Create folder to store later results (te: time evolution)
    os.makedirs('data_post_te', mode=0o777, exist_ok=True)

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
    
    # Initialize the array for alias of mean streamwise velocity profile averaged with different flow realizations
    # at a certain time unit
    mean_u = np.zeros(ny)

    # Initialize instantaneous mean statistics array, not averaged with different flow realizations
    # we have 9 different statistics (mean, var, Reynolds stress) (function of y)
    mean_stats = np.zeros((ny, 9))
    
    # Initialize mean statistics array, function of y and time
    mean_stats_realiz = np.zeros((ny, 9, nsavings))

    # Initialize arrays for TTBL thickness parameters
    delta_99 = np.zeros(nsavings)   # BL thickness delta_99
    disp_t   = np.zeros(nsavings)   # displacement thickness, delta*
    mom_t    = np.zeros(nsavings)   # momentum     thickness, theta

    """
    Do loop over all savings of 'mean_stats_runtime' files. 
    The ts goes from the first saving (IC, ts = 1) to the last one, ilast, with increment ioutput_cf.
    Index for the cycle is ti: time index.
    """
    
    #!--- Save averaged mean statistics ---!
    print()
    print(">>> Saving 'mean_stats_realiz' in /data_post_te.")
    print()
        
    for ti in range(0, nsavings, 1):
      
        #!--- Calculate ts to open 'mean_stats_runtime-ts' file (ts_iter: time-step of the iterations) ---!
        
        # ts of the initial condition is ts = 1
        if ti == 0:

            ts_iter = 1

        # All the other ts are multiples of 'output_cf'
        else:

            ts_iter = ti*ioutput_cf
        
        # Do loop over different realizations, from 1 to nr
        for i in range(1, nr+1, 1):
               
            # Read of mean statistics calculated runtime data from 'data/mean_stats_runtime' folder
            mean_stats = np.loadtxt(f'data_r{i:01d}/mean_stats_runtime/mean_stats_runtime-ts{ts_iter:07d}.txt', skiprows=12, delimiter=',', dtype=np.float64)
                        
            # Summing mean statistics array with different realizations into the overall array for time-evolution
            mean_stats_realiz[:,:,ti] = mean_stats_realiz[:,:,ti] + mean_stats[:,:] / nr
            
        #!--- Finalize 2nd order statistics ---!
        
        # Variances
        mean_stats_realiz[:,3,ti] = mean_stats_realiz[:,3,ti] - mean_stats_realiz[:,0,ti]**2  # streamwise  velocity variance
        mean_stats_realiz[:,4,ti] = mean_stats_realiz[:,4,ti] - mean_stats_realiz[:,1,ti]**2  # wall-normal velocity variance
        mean_stats_realiz[:,5,ti] = mean_stats_realiz[:,5,ti] - mean_stats_realiz[:,2,ti]**2  # spanwise    velocity variance    
        
        # Reynolds stress
        mean_stats_realiz[:,6,ti] = mean_stats_realiz[:,6,ti] - mean_stats_realiz[:,0,ti]*mean_stats_realiz[:,1,ti]  # Reynolds stress <u'v'>
        mean_stats_realiz[:,7,ti] = mean_stats_realiz[:,7,ti] - mean_stats_realiz[:,0,ti]*mean_stats_realiz[:,2,ti]  # Reynolds stress <u'w'>
        mean_stats_realiz[:,8,ti] = mean_stats_realiz[:,8,ti] - mean_stats_realiz[:,1,ti]*mean_stats_realiz[:,2,ti]  # Reynolds stress <v'w'>

        # Create the file and write; we are adding at each file the shear velocities coming from the main function 
        with open(f'data_post_te/mean_stats_realiz-ts{ts_iter:07d}.txt', 'w') as f:
            f.write(f'Mean statistics at ts={ts_iter}.\n')        
            f.write('\n')
            f.write(f"{'sh_vel_x':>{pp.c_w}} = {sh_velx[ti]:{pp.fs6}}\n")
            f.write(f"{'sh_vel_tot':>{pp.c_w}} = {sh_veltot[ti]:{pp.fs6}}\n")
            f.write('\n') 
            f.write(f"{'mean[u]':>{pp.c_w}}, "    +
                    f"{'mean[v]':>{pp.c_w}}, "    +
                    f"{'mean[w]':>{pp.c_w}}, "    +
                    f"{'var[u]':>{pp.c_w}}, "     +
                    f"{'var[v]':>{pp.c_w}}, "     +
                    f"{'var[w]':>{pp.c_w}}, "     +
                    f"{'mean[uv]':>{pp.c_w}}, "   +
                    f"{'mean[uw]':>{pp.c_w}}, "   +
                    f"{'mean[vw]':>{pp.c_w}}\n"   )
                
            for j in range(0, ny):
                f.write(f"{mean_stats_realiz[j,0,ti]:{pp.fs6}}, " +
                        f"{mean_stats_realiz[j,1,ti]:{pp.fs6}}, " +
                        f"{mean_stats_realiz[j,2,ti]:{pp.fs6}}, " +
                        f"{mean_stats_realiz[j,3,ti]:{pp.fs6}}, " +
                        f"{mean_stats_realiz[j,4,ti]:{pp.fs6}}, " +
                        f"{mean_stats_realiz[j,5,ti]:{pp.fs6}}, " +
                        f"{mean_stats_realiz[j,6,ti]:{pp.fs6}}, " +
                        f"{mean_stats_realiz[j,7,ti]:{pp.fs6}}, " +
                        f"{mean_stats_realiz[j,8,ti]:{pp.fs6}}\n" )

        #!--- Calculation of thickness parameters ---!
        
        # Take alias for mean streamwise velocity profile
        mean_u[:] = mean_stats_realiz[:,0,ti]
               
        # Calculate BL thickness delta_99 for a TTBL and its related index
        (bl_thick, bl_thick_j) = calculate_ttbl_delta_99(mean_u, yp)
        
        # Store result in 'delta_99' array
        delta_99[ti] = bl_thick
        
        # Calculate the displacement thickness delta*
        int1 = mean_u/uwall  # 'integrand 1' 

        # Interpolation at the 6th order of accuracy with a spline of 5th order, 
        # that passes through all data points
        spl = InterpolatedUnivariateSpline(yp, int1, k=5)
        disp_t[ti] = spl.integral(y0, yn)

        # Calculate the momentum thickness theta
        int2 = int1 - int1**2  # 'integrand 2' 

        # Interpolation at the 6th order of accuracy with a spline of 5th order,
        # that passes through all data points
        spl = InterpolatedUnivariateSpline(yp, int2, k=5)
        mom_t[ti] = spl.integral(y0, yn)

    # Return to main program
    return (delta_99, disp_t, mom_t)



