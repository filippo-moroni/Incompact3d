
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

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '../../4-common', 'python_common'))
sys.path.append(config_path)

# Import the plotting_params module
import plot_params as pp

# Import function to read 'input.i3d'
from read_files import read_input_files

# Import function to setup flow parameters and mesh 
from set_flow_parameters import set_flow_parameters

# Import function to calculate TTBL thickness parameters
from ttbl_subs import calculate_ttbl_thick_params

# Import shutil to delete folders with contents inside and subfolders
import shutil

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this subroutine we perform:
!               - averaging between different flow realizations of runtime 
!                 mean statistics files 'mean_stats_runtime' and saving of 
!                 the related shear velocities in the average files 
!                 'mean_stats_realiz-ts'; 
!               - 6th order accurate calculations of integrals of TTBL 
!                 thickness parameters (delta*, theta) and delta_99 using 
!                 an external subroutine applied on the flow field averaged 
!                 with different flow realizations.                              
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def average_runtime_mean_stats(sh_vel_tot, sh_vel_x, mg_phi_w, nsavings, time_window_index, den):

    """
    Inputs: 
    
     - sh_vel_tot        : total shear velocity;
     - sh_vel_x          : streamwise shear velocity;
     - mg_phi_w          : mean scalar gradient at the wall;
     - nsavings          : number of savings in time (defined in the main function);
     - time_window_index : number of savings used to time-window average the temporal snapshots;
                           this index refers to half of the window excluding the central snapshot
                           (e.g. if == 1: a window of 3 snapshots);
     - den               : denominator of the divisions (number of snapshots in the time window average times
                           number of flow realizations).      
    
    Outputs:
    
     - mean statistics saved runtime averaged with different flow realizations;
     - mean scalar statistics saved runtime averaged with different flow realizations; 
     - TTBL thickness parameters to check the flow evolution.
    """
    
    # Delete the folder of time-evolution results to be sure we are looking at new results
    shutil.rmtree('data_post_te', ignore_errors=not os.path.exists('data_post_te'))
    
    # Create folder to store later results (te: time evolution)
    os.makedirs('data_post_te',          mode=0o777, exist_ok=True)
    os.makedirs('data_post_te/velocity', mode=0o777, exist_ok=True)
    os.makedirs('data_post_te/scalar',   mode=0o777, exist_ok=True)

    #!--- Reading of files section and setup of flow parameters ---!

    # Read useful flow parameters from 'input.i3d' and 'post.prm' files
    (itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, itimescheme, ioutput, ioutput_cf, iswitch_wo,  
     add_string, file1, filen, icrfile, nr, post_mean, post_grad, post_corz, post_tke_eq
    ) = read_input_files('input.i3d','post.prm')

    #!--- Parameters and mesh ---!
    (nu, y, uwall, twd, phiwall) = set_flow_parameters(re)
                          
    #!--------------------------------------------------------------------------------------!
                                    
    """
    !---------------------------------------------------------------------------------------------------------------------!
     Arrays for time-evolution of (instantaneous) mean statistics of a specific flow realization.
     These statistics are averaged only in homogeneous directions (x and z).
     
      - ny       : number of points in y-direction; 
      - nsavings : number of savings in time due to the Incompact3d 'modified' solver 'print_cf' (see 'extra_tools.f90').
                   This is an external input of the subroutine, calculated from 'cf_history.txt'. 
    !---------------------------------------------------------------------------------------------------------------------!
    """
    
    mean_u  = np.zeros((ny, nsavings))     # mean streamwise  velocity
    mean_v  = np.zeros((ny, nsavings))     # mean wall-normal velocity
    mean_w  = np.zeros((ny, nsavings))     # mean spanwise    velocity
    var_u   = np.zeros((ny, nsavings))     # streamwise  variance
    var_v   = np.zeros((ny, nsavings))     # wall-normal variance
    var_w   = np.zeros((ny, nsavings))     # spanwise    variance
    mean_uv = np.zeros((ny, nsavings))     # Reynolds stress <u'v'>
    mean_uw = np.zeros((ny, nsavings))     # Reynolds stress <u'w'>
    mean_vw = np.zeros((ny, nsavings))     # Reynolds stress <v'w'>
    
    """
    !---------------------------------------------------------------------------------------------------------------------!
     Arrays for time-evolution of mean statistics, averaged with different flow realizations and time-window average 
     if requested (r: averaged on flow realizations).
     
      - ny           : number of points in y-direction; 
      - nsavings_red : number of savings in time; if time-window average is requested, the first and the last savings of 
                       the entire time evolution are not stored, since the time-window average is centered on a generic 
                       time step ts . 
    !---------------------------------------------------------------------------------------------------------------------!
    """
    
    # Number of elements in the reduced arrays (number of savings reduced)
    nsavings_red = nsavings - 2*time_window_index
    
    mean_u_r  = np.zeros((ny, nsavings_red))  # mean streamwise  velocity, flow realizations and time-window averaged
    mean_v_r  = np.zeros((ny, nsavings_red))  # mean wall-normal velocity, flow realizations and time-window averaged
    mean_w_r  = np.zeros((ny, nsavings_red))  # mean spanwise    velocity, flow realizations and time-window averaged
    var_u_r   = np.zeros((ny, nsavings_red))  # streamwise  variance,      flow realizations and time-window averaged
    var_v_r   = np.zeros((ny, nsavings_red))  # wall-normal variance,      flow realizations and time-window averaged
    var_w_r   = np.zeros((ny, nsavings_red))  # spanwise    variance,      flow realizations and time-window averaged
    mean_uv_r = np.zeros((ny, nsavings_red))  # Reynolds stress <u'v'>,    flow realizations and time-window averaged
    mean_uw_r = np.zeros((ny, nsavings_red))  # Reynolds stress <u'w'>,    flow realizations and time-window averaged
    mean_vw_r = np.zeros((ny, nsavings_red))  # Reynolds stress <v'w'>,    flow realizations and time-window averaged
        
    # Initialize scalar arrays if present
    if numscalar == 1:
            
        mean_phi  = np.zeros((ny, nsavings))  # mean scalar field
        var_phi   = np.zeros((ny, nsavings))  # scalar variance
        mean_uphi = np.zeros((ny, nsavings))  # mixed fluctuations <u'phi'>
        mean_vphi = np.zeros((ny, nsavings))  # mixed fluctuations <v'phi'>
        mean_wphi = np.zeros((ny, nsavings))  # mixed fluctuations <w'phi'>
        
        mean_phi_r  = np.zeros((ny, nsavings_red))  # mean scalar field,           flow realizations and time-window averaged
        var_phi_r   = np.zeros((ny, nsavings_red))  # scalar variance,             flow realizations and time-window averaged
        mean_uphi_r = np.zeros((ny, nsavings_red))  # mixed fluctuations <u'phi'>, flow realizations and time-window averaged
        mean_vphi_r = np.zeros((ny, nsavings_red))  # mixed fluctuations <v'phi'>, flow realizations and time-window averaged
        mean_wphi_r = np.zeros((ny, nsavings_red))  # mixed fluctuations <w'phi'>, flow realizations and time-window averaged
        
        
    # Initialize arrays for TTBL thickness parameters
    delta_99   = np.zeros(nsavings_red)   # BL thickness delta_99
    disp_t     = np.zeros(nsavings_red)   # displacement thickness, delta*
    mom_t      = np.zeros(nsavings_red)   # momentum     thickness, theta
    
    # Initialize maximum Delta y+ at the BL edge for grid resolutions monitoring
    max_delta_yd_plus = 0.0

    """
    !---------------------------------------------------------------------------------------------------------------------!
     Do loop over all savings of 'mean_stats_runtime' files. 
     The ts goes from the first saving (IC, ts = 1) to the last one, ilast, with increment ioutput_cf.
     Index for the cycle is ti: time index.
    !---------------------------------------------------------------------------------------------------------------------!
    
    """
    
    #!--- Average mean statistics with different flow realizations ---!
    print()
    print(">>> Averaging 'mean_stats_runtime' files with different flow realizations.")
    print()
    
    # Do loop over different realizations, from 1 to nr
    for i in range(1, nr+1, 1):
    
        # Do loop from 0 to number of savings (ti: time index, that represents the different savings in time)    
        for ti in range(0, nsavings, 1):
      
            #!--- Calculate ts to open 'mean_stats_runtime-ts' file (ts_iter: time-step of the iterations) ---!
        
            # ts of the initial condition is ts = 1
            if ti == 0: ts_iter = 1

            # All the other ts are multiples of 'output_cf'
            else: ts_iter = ti*ioutput_cf
                
            # Read of mean statistics calculated runtime data from 'data/mean_stats_runtime/velocity' folder
            mean_stats = np.loadtxt(f'data_r{i:01d}/mean_stats_runtime/velocity/mean_stats_runtime-ts{ts_iter:07d}.txt', skiprows=12, delimiter=',', dtype=np.float64)
            
            mean_u [:,ti] = mean_stats[:,0]
            mean_v [:,ti] = mean_stats[:,1]
            mean_w [:,ti] = mean_stats[:,2]
            var_u  [:,ti] = mean_stats[:,3]
            var_v  [:,ti] = mean_stats[:,4]           
            var_w  [:,ti] = mean_stats[:,5]
            mean_uv[:,ti] = mean_stats[:,6]
            mean_uw[:,ti] = mean_stats[:,7]
            mean_vw[:,ti] = mean_stats[:,8]            
            
            # Read of mean scalar statistics calculated runtime data from 'data/mean_stats_runtime/scalar' folder
            if numscalar == 1:
                       
                mean_stats_scalar = np.loadtxt(f'data_r{i:01d}/mean_stats_runtime/scalar/mean_stats_scalar_runtime-ts{ts_iter:07d}.txt', skiprows=9, delimiter=',', dtype=np.float64)
            
                mean_phi [:,ti] = mean_stats_scalar[:,0]
                var_phi  [:,ti] = mean_stats_scalar[:,1]
                mean_uphi[:,ti] = mean_stats_scalar[:,2]
                mean_vphi[:,ti] = mean_stats_scalar[:,3]               
                mean_wphi[:,ti] = mean_stats_scalar[:,4]
                
        # Take alias for time-window index
        twi = time_window_index        
        
        # Cycle to sum with time-window average
        for ti in range(0+twi, nsavings-twi, 1):
            
            # Time-window average cycle
            for i in range(-twi, twi+1, 1):
                        
                # Summing mean statistics arrays with different realizations into the overall arrays for time-evolution
                mean_u_r [:,ti-twi] = mean_u_r [:,ti-twi] + mean_u [:,ti+i] / den
                mean_v_r [:,ti-twi] = mean_v_r [:,ti-twi] + mean_v [:,ti+i] / den
                mean_w_r [:,ti-twi] = mean_w_r [:,ti-twi] + mean_w [:,ti+i] / den
                var_u_r  [:,ti-twi] = var_u_r  [:,ti-twi] + var_u  [:,ti+i] / den
                var_v_r  [:,ti-twi] = var_v_r  [:,ti-twi] + var_v  [:,ti+i] / den                
                var_w_r  [:,ti-twi] = var_w_r  [:,ti-twi] + var_w  [:,ti+i] / den
                mean_uv_r[:,ti-twi] = mean_uv_r[:,ti-twi] + mean_uv[:,ti+i] / den           
                mean_uw_r[:,ti-twi] = mean_uw_r[:,ti-twi] + mean_uw[:,ti+i] / den            
                mean_vw_r[:,ti-twi] = mean_vw_r[:,ti-twi] + mean_vw[:,ti+i] / den            
            
                # Summing mean scalar statistics arrays with different realizations into the overall arrays for time-evolution
                if numscalar == 1:
                
                    mean_phi_r [:,ti-twi] = mean_phi_r [:,ti-twi] + mean_phi [:,ti+i] / den
                    var_phi_r  [:,ti-twi] = var_phi_r  [:,ti-twi] + var_phi  [:,ti+i] / den                
                    mean_uphi_r[:,ti-twi] = mean_uphi_r[:,ti-twi] + mean_uphi[:,ti+i] / den
                    mean_vphi_r[:,ti-twi] = mean_vphi_r[:,ti-twi] + mean_vphi[:,ti+i] / den
                    mean_wphi_r[:,ti-twi] = mean_wphi_r[:,ti-twi] + mean_wphi[:,ti+i] / den    
    
    #!--- Save averaged mean statistics ---!
    print(">>> Finalizing and saving 'mean_stats_realiz-ts' files in /data_post_te.")
    print()
                    
    # Do loop from 0 to number of savings (ti: time index, that represents the different savings in time) 
    for ti in range(0, nsavings_red, 1):
    
        # Time-step of central saving
        ts_iter = ti*ioutput_cf    
                               
        #!--- Finalize 2nd order statistics ---!
        
        # Variances
        var_u_r[:,ti] = var_u_r[:,ti] - mean_u_r[:,ti]**2  # streamwise  velocity variance
        var_v_r[:,ti] = var_v_r[:,ti] - mean_v_r[:,ti]**2  # wall-normal velocity variance
        var_w_r[:,ti] = var_w_r[:,ti] - mean_w_r[:,ti]**2  # spanwise    velocity variance 
        
        # Reynolds stress
        mean_uv_r[:,ti] = mean_uv_r[:,ti] - mean_u_r[:,ti]*mean_v_r[:,ti]  # Reynolds stress <u'v'>
        mean_uw_r[:,ti] = mean_uw_r[:,ti] - mean_u_r[:,ti]*mean_w_r[:,ti]  # Reynolds stress <u'w'>
        mean_vw_r[:,ti] = mean_vw_r[:,ti] - mean_v_r[:,ti]*mean_w_r[:,ti]  # Reynolds stress <v'w'>
                                
        # Scalar field
        if numscalar == 1:
        
            # Variance
            var_phi_r[:,ti] = var_phi_r[:,ti] - mean_phi_r[:,ti]**2  # scalar field variance
                        
            # Mixed fluctuations
            mean_uphi_r[:,ti] = mean_uphi_r[:,ti] - mean_u_r[:,ti]*mean_phi_r[:,ti]  # mixed fluctuations <u'phi'>
            mean_vphi_r[:,ti] = mean_vphi_r[:,ti] - mean_v_r[:,ti]*mean_phi_r[:,ti]  # mixed fluctuations <v'phi'>            
            mean_wphi_r[:,ti] = mean_wphi_r[:,ti] - mean_w_r[:,ti]*mean_phi_r[:,ti]  # mixed fluctuations <w'phi'>
            
        # Create the file and write; we are adding at each file the shear velocities coming from the main function 
        with open(f'data_post_te/velocity/mean_stats_realiz-ts{ts_iter:07d}.txt', 'w') as f:
            f.write(f'Mean statistics at ts={ts_iter}.\n')        
            f.write('\n')
            f.write(f"{'sh_vel_x':>{pp.c_w}} = {sh_vel_x[ti]:{pp.fs6}}\n")
            f.write(f"{'sh_vel_tot':>{pp.c_w}} = {sh_vel_tot[ti]:{pp.fs6}}\n")
            f.write('\n') 
            f.write(f"{'mean[u]':>{pp.c_w}}, "  +
                    f"{'mean[v]':>{pp.c_w}}, "  +
                    f"{'mean[w]':>{pp.c_w}}, "  +
                    f"{'var[u]':>{pp.c_w}}, "   +
                    f"{'var[v]':>{pp.c_w}}, "   +
                    f"{'var[w]':>{pp.c_w}}, "   +
                    f"{'mean[uv]':>{pp.c_w}}, " +
                    f"{'mean[uw]':>{pp.c_w}}, " +
                    f"{'mean[vw]':>{pp.c_w}}\n" )
                
            for j in range(0, ny):
                f.write(f"{mean_u_r[j,ti]:{pp.fs6}}, "  +
                        f"{mean_v_r[j,ti]:{pp.fs6}}, "  +
                        f"{mean_w_r[j,ti]:{pp.fs6}}, "  +
                        f"{var_u_r[j,ti]:{pp.fs6}}, "   +
                        f"{var_v_r[j,ti]:{pp.fs6}}, "   +
                        f"{var_w_r[j,ti]:{pp.fs6}}, "   +
                        f"{mean_uv_r[j,ti]:{pp.fs6}}, " +
                        f"{mean_uw_r[j,ti]:{pp.fs6}}, " +
                        f"{mean_vw_r[j,ti]:{pp.fs6}}\n" )
        
        # Save runtime statistics averaged on flow realizations for scalar field
        if numscalar == 1:
                
            # Create the file and write; we are adding at each file the mean scalar gradient at the wall coming from the main function 
            with open(f'data_post_te/scalar/mean_stats_scalar_realiz-ts{ts_iter:07d}.txt', 'w') as f:
                f.write(f'Mean scalar statistics at ts={ts_iter}.\n')        
                f.write('\n')
                f.write(f"{'(dPhi/dy)w':>{pp.c_w}} = {mg_phi_w[ti]:{pp.fs6}}\n")
                f.write('\n') 
                f.write(f"{'mean[phi]':>{pp.c_w}}, " +
                        f"{'var[phi]':>{pp.c_w}}, "  +
                        f"{'<u phi>':>{pp.c_w}}, "   +
                        f"{'<v phi>':>{pp.c_w}}, "   +
                        f"{'<w phi>':>{pp.c_w}}\n"   )
                
                for j in range(0, ny):
                    f.write(f"{mean_phi_r[j,ti]:{pp.fs6}}, "  +
                            f"{var_phi_r[j,ti]:{pp.fs6}}, "   +
                            f"{mean_uphi_r[j,ti]:{pp.fs6}}, " +
                            f"{mean_vphi_r[j,ti]:{pp.fs6}}, " +
                            f"{mean_wphi_r[j,ti]:{pp.fs6}}\n" )
      
      
        #!--- Calculation of thickness parameters ---!
                
        # Call of external subroutine for the calculation of TTBL thickness parameters
        (delta_99[ti], delta_99_j, disp_t[ti], mom_t[ti]) = calculate_ttbl_thick_params(mean_u_r[:,ti],y,uwall)
        
        # Delta y+ at the BL edge
        delta_yd_plus = (y[delta_99_j] - y[delta_99_j-1]) * sh_vel_tot[ti] / nu
        
        # Maximum Delta y+ at the BL edge for grid resolutions monitoring
        if max_delta_yd_plus < delta_yd_plus:
            max_delta_yd_plus = delta_yd_plus
               
    # Return to main program with TTBL thickness parameters time evolution and maximum delta y+ at the BL interface
    return (delta_99, disp_t, mom_t, max_delta_yd_plus)
    



