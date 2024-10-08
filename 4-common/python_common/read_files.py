
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: In this file, we store Python functions to read .txt files 
!              for post-processing and plots:                             
!               - read_input_files:        to read 'input.i3d' and 'post.prm';    
!               - read_data:               to read statistics data;               
!               - read_ref_data:           to read reference data;
!               - read_ref_data_temp_evol: to read TTBL temporal evolution 
!                                          quantities.  
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Common libraries
import numpy as np
import os

# Get environment variable for Desktop path
base_path = os.getenv('desktop')

if base_path is None:
    raise EnvironmentError("Please set the Desktop environment variable.")

# Directory where reference data are stored
dirname = os.path.join(base_path, 'Incompact3d/3-post_processing/reference_data')

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Function to read Incompact3d files: 'input.i3d' and 'post.prm'.   
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def read_input_files(filename1,filename2):
      
    # Opening of 'input.i3d' file
    with open(filename1, 'r') as file:
        
        # Read all lines into a list
        lines = file.readlines()
    
    # Extract: itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, 
    # Re, dt, ifirst, ilast, numscalar, ioutput, ioutput_cf, iswitch_wo 
        
    # As always, index is 1 less of the line number (Python convention)
    itype       = lines[7]  
    nx          = lines[14]
    ny          = lines[15]
    nz          = lines[16]
    istret      = lines[17]
    beta        = lines[18]
    Lx          = lines[21]
    Ly          = lines[22]
    Lz          = lines[23]
    re          = lines[26]
    dt          = lines[29]
    ifirst      = lines[30]
    ilast       = lines[31]
    numscalar   = lines[35]
    itimescheme = lines[61]
    ioutput     = lines[73]
    ioutput_cf  = lines[74]
    iswitch_wo  = lines[90]
    
    # Removing characters in front of the extracted strings and the comments:
    # 1) split: the string is split when the specified character is encountered; 
    # 2) we select the portion of string with index inside square brackets;
    # 3) strip: removes leading or trailing whitespaces from the string. 
    
    itype       = itype.split('=')[-1].strip()
    
    nx          = nx.split('!')[0]
    nx          = nx.split('=')[-1].strip()
        
    ny          = ny.split('!')[0]
    ny          = ny.split('=')[-1].strip()
    
    nz          = nz.split('!')[0]
    nz          = nz.split('=')[-1].strip()
        
    istret      = istret.split('!')[0]
    istret      = istret.split('=')[-1].strip()
        
    beta        = beta.split('!')[0]
    beta        = beta.split('=')[-1].strip()
    
    Lx          = Lx.split('!')[0]
    Lx          = Lx.split('=')[-1].strip()
    
    Ly          = Ly.split('!')[0]
    Ly          = Ly.split('=')[-1].strip()
    
    Lz          = Lz.split('!')[0]
    Lz          = Lz.split('=')[-1].strip()
    
    re          = re.split('!')[0]
    re          = re.split('=')[-1].strip()
        
    dt          = dt.split('!')[0]
    dt          = dt.split('=')[-1].strip()
        
    ifirst      = ifirst.split('!')[0]
    ifirst      = ifirst.split('=')[-1].strip()
        
    ilast       = ilast.split('!')[0]
    ilast       = ilast.split('=')[-1].strip()
        
    numscalar   = numscalar.split('!')[0]
    numscalar   = numscalar.split('=')[-1].strip()
    
    itimescheme = itimescheme.split('!')[0]
    itimescheme = itimescheme.split('=')[-1].strip()
        
    ioutput     = ioutput.split('!')[0]
    ioutput     = ioutput.split('=')[-1].strip()
        
    ioutput_cf  = ioutput_cf.split('!')[0]
    ioutput_cf  = ioutput_cf.split('=')[-1].strip()
    
    iswitch_wo  = iswitch_wo.split('!')[0]
    iswitch_wo  = iswitch_wo.split('=')[-1].strip()
    
    # Convert to needed variable type (integer, float, etc.)
    itype       = int(itype)
    nx          = int(nx)
    ny          = int(ny)
    nz          = int(nz)
    istret      = int(istret)
    beta        = np.float64(beta)
    Lx          = np.float64(Lx)
    Ly          = np.float64(Ly)
    Lz          = np.float64(Lz)
    re          = np.float64(re)
    dt          = np.float64(dt)
    ifirst      = int(ifirst)
    ilast       = int(ilast)
    numscalar   = int(numscalar)
    itimescheme = int(itimescheme)
    ioutput     = int(ioutput)
    ioutput_cf  = int(ioutput_cf)
    iswitch_wo  = int(iswitch_wo)
    
    # Opening of 'post.prm' file
    with open(filename2, 'r') as file:
        
        # Read all lines into a list
        lines = file.readlines()
   
    # Extract needed lines  
    add_string  = lines[3]    # Flow case name
    file1       = lines[5]    # First snapshot index
    filen       = lines[6]    # Final snapshot index
    icrfile     = lines[7]    # File increment
    nr          = lines[8]    # Number of flow realizations
        
    post_mean   = lines[12]   # Compute mean statistics
    post_vort   = lines[13]   # Compute mean vorticity and mean gradient
    post_diss   = lines[14]   # Compute mean total dissipation rate    
    post_corz   = lines[15]   # Compute correlation functions along z (a previous run with post_mean = 1 must be performed)
    post_tke_eq = lines[16]   # Compute fluctuating terms of TKE equations (a previous run with post_mean = 1 must be performed)
        
    # Extract the needed variables
    add_string  = add_string.split('!')[0]
    add_string  = add_string.rstrip()
        
    file1       =   file1.split('#')[0].strip()
    filen       =   filen.split('#')[0].strip()
    icrfile     = icrfile.split('#')[0].strip()
    nr          =      nr.split('#')[0].strip()
        
    post_mean   =   post_mean.split('#')[0].strip()
    post_vort   =   post_vort.split('#')[0].strip()
    post_diss   =   post_diss.split('#')[0].strip()
    post_corz   =   post_corz.split('#')[0].strip()
    post_tke_eq = post_tke_eq.split('#')[0].strip()
        
    # Convert to needed variable type (integer or boolean)
    file1       = int(file1)
    filen       = int(filen)
    icrfile     = int(icrfile)
    nr          = int(nr)
        
    post_mean   = bool(int(post_mean))
    post_vort   = bool(int(post_vort))
    post_diss   = bool(int(post_diss))
    post_corz   = bool(int(post_corz))
    post_tke_eq = bool(int(post_tke_eq))
    
    # Halve the points in y direction for a channel
    if itype == 3:
        ny = (ny - 1) // 2 + 1
        
    # Return to main program with extracted parameters
    return (
            # From 'input.i3d'
            itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, itimescheme, ioutput, ioutput_cf, iswitch_wo,
            
            # From 'post.prm' 
            add_string, file1, filen, icrfile, nr, post_mean, post_vort, post_diss, post_corz, post_tke_eq
           ) 

#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Python function to read statistics data, obtained from
!              'post_incompact3d'. For TTBL simulations only, we can plot
!              statistics from 'cf_monitoring'. 
!              For data from 'post_incompact3d', shear velocities are also 
!              calculated from mean gradients. On the other hand, 
!              'cf_monitoring' data contains already directly the shear 
!              velocities.
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""
    
def read_data(itype, numscalar, post_mean, post_vort, post_diss, post_corz, post_tke_eq, ny, nz, nu, snap_numb=None):

    #!--- Initialize variables ---! 
    
    # Snapshot number
    snap_numb   = None
    
    # Mean statistics
    mean_u      = 0.0
    mean_w      = 0.0
    var_u       = 0.0
    var_v       = 0.0
    var_w       = 0.0
    mean_uv     = 0.0
    
    # Mean vorticity and mean gradients
    vort_x      = 0.0   # streamwise mean vorticity
    vort_y      = 0.0   # wall-normal mean vorticity
    vort_z      = 0.0   # spanwise mean vorticity
    mg_x        = 0.0   # streamwise mean gradient
    mg_z        = 0.0   # spanwise mean gradient
    mg_phi      = 0.0   # scalar mean gradient
    
    # Total dissipation
    eps         = 0.0
    
    # Spanwise correlation functions
    Ruuz = np.zeros((ny,nz), dtype=np.float64, order='F')
    Rvvz = np.zeros((ny,nz), dtype=np.float64, order='F')
    Rwwz = np.zeros((ny,nz), dtype=np.float64, order='F')
    Ruvz = np.zeros((ny,nz), dtype=np.float64, order='F')
    Rssz = np.zeros((ny,nz), dtype=np.float64, order='F')
    
    # Turbulent Kinetic Energy (TKE) budget
    tke_turbt  = 0.0   
    tke_presst = 0.0     
    tke_difft  = 0.0
    tke_prod   = 0.0
    tke_pseps  = 0.0
    
    # Shear velocities
    sh_vel_x   = 0.0  
    sh_vel_tot = 0.0
    
    # Scalar field friction quantities
    phi_tau = 0.0
    
    # Time-step
    ts = None
    
    # Switcher to select the source of statistics for a TTBL (0: from 'post_incompact3d', 1: from 'cf_monitoring')
    i_switch_plot = None
    
    print()

    # Channel
    if itype == 3:
    
        print(">>> Plotting of statistics for a channel.")
        
        # Channel statistics are only from 'post_incompact3d' at the moment
        i_switch_plot = False

        # Reading of mean statistics
        if post_mean:        
            M = np.loadtxt('data_post/mean_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
            mean_u  =   M[:,0]
            mean_w  =   M[:,2]   
            var_u   =   M[:,3]
            var_v   =   M[:,4]
            var_w   =   M[:,5]
            mean_uv = - M[:,12]  # Change sign for Reynolds stress <u'v'> for a Channel
    
        """
        Reading of vorticity components and mean gradients
        (always performed since we need the mean gradients to calculate
        total shear velocity).
        
        """
        M = np.loadtxt('data_post/vort_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
        vort_x = M[:,0]
        vort_y = M[:,1]
        vort_z = M[:,2]
        mg_x   = M[:,3]
        mg_z   = M[:,4]
            
        # Reading of the mean total dissipation
        if post_diss:
            eps = np.loadtxt('data_post/diss_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
        # Reading of correlations
        if post_corz:
            Ruuz = np.loadtxt('data_post/Ruuz.txt', skiprows=0, delimiter=None, dtype=np.float64)
            Rvvz = np.loadtxt('data_post/Rvvz.txt', skiprows=0, delimiter=None, dtype=np.float64)
            Rwwz = np.loadtxt('data_post/Rwwz.txt', skiprows=0, delimiter=None, dtype=np.float64)
            Ruvz = np.loadtxt('data_post/Ruvz.txt', skiprows=0, delimiter=None, dtype=np.float64)
      
            # Read scalar field correlations
            if numscalar == 1:               
                Rssz = np.loadtxt('data_post/Rssz.txt', skiprows=0, delimiter=None, dtype=np.float64)
        
        # Read of TKE 
        if post_tke_eq:
            M = np.loadtxt('data_post/tke_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
            tke_turbt  = M[:,0]   
            tke_presst = M[:,1]     
            tke_difft  = M[:,2]
            tke_prod   = M[:,3]
            tke_pseps  = M[:,4]
        
    # TTBL
    elif itype == 13:

        print(">>> Plotting of statistics for a TTBL.")
        print()

        """
        Asking to the user if he wants to plot statistics calculated from:
         - 'post_incompact3d', that employs full snapshots;
         - 'cf_monitoring', that employs mean statistics calculated runtime;
            in this case, only mean statistics can be plotted. 
        """
        print(">>> Specify the type of results to plot: ")
        print("    - 0: from 'post_incompact3d' ")  
        print("    - 1: from 'cf_monitoring' ")          
        print()
        
        # Switcher for plotting different sources of statistics
        i_switch_plot = bool(int(input()))

        # Plotting statistics from 'post_incompact3d'
        if i_switch_plot == False:
        
            # Asking to the user the specific snapshot to show
            snap_numb = input(">>> Enter the snapshot number to show: ")
    
            # Pad with zeros to match snapshots' naming  
            snap_numb = snap_numb.zfill(4)
                     
            """
            Reading of mean statistics. For a TTBL it is always done since 
            we need to calculate the boundary layer thickness delta_99 to display
            friction Reynolds number.
            """
            M = np.loadtxt(f'data_post/mean_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
            mean_u  = M[:,0]
            mean_w  = M[:,2]   
            var_u   = M[:,3]
            var_v   = M[:,4]
            var_w   = M[:,5]
            mean_uv = M[:,12]
    
            """
            Reading of vorticity components and mean gradients
            (always performed since we need the mean gradients to calculate
            total shear velocity).
            """
            M = np.loadtxt(f'data_post/vort_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
            vort_x = M[:,0]
            vort_y = M[:,1]
            vort_z = M[:,2]
            mg_x   = M[:,3]
            mg_z   = M[:,4]
            mg_phi = M[:,5]
                
            # Reading of the mean total dissipation
            if post_diss:
                eps = np.loadtxt(f'data_post/diss_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
        
            # Reading of correlations
            if post_corz:
                Ruuz = np.loadtxt(f'data_post/Ruuz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)
                Rvvz = np.loadtxt(f'data_post/Rvvz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)
                Rwwz = np.loadtxt(f'data_post/Rwwz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)
                Ruvz = np.loadtxt(f'data_post/Ruvz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)
        
                # Read scalar field correlations
                if numscalar == 1:
                    Rssz = np.loadtxt(f'data_post/Rssz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)
        
            # Read of TKE
            if post_tke_eq: 
                M = np.loadtxt(f'data_post/tke_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
                tke_turbt  = M[:,0]   
                tke_presst = M[:,1]     
                tke_difft  = M[:,2]
                tke_prod   = M[:,3]
                tke_pseps  = M[:,4]

        # Plotting statistics from 'cf_monitoring'
        elif i_switch_plot == True:
        
            # Asking to the user the time-step to show
            ts = input(">>> Enter the time-step to show: ")
            
            # Pad with zeros to match 'mean_stats_realiz-ts' naming  
            ts = ts.zfill(7)
           
            # Read mean stats averaged with different flow realizations.
            # These are time-dependent statistics.
            M = np.loadtxt(f'data_post_te/velocity/mean_stats_realiz-ts{ts}.txt', skiprows=6, delimiter=',', dtype=np.float64)
                        
            mean_u  = M[:,0]
            mean_w  = M[:,2]   
            var_u   = M[:,3]
            var_v   = M[:,4]
            var_w   = M[:,5]
            mean_uv = M[:,6]
            
            # Opening of 'mean_stats_realiz-ts' file
            with open(f'data_post_te/velocity/mean_stats_realiz-ts{ts}.txt', 'r') as file:
    
                # Read all lines into a list
                lines = file.readlines()
            
            # Shear velocities reading: as always, index is 1 less of the line number (Python convention)
            sh_vel_x   = lines[2]   # streamwise shear velocity (based on streamwise mean gradient)  
            sh_vel_tot = lines[3]   # total shear velocity (based on total mean gradient) 
            
            # Removing characters in front of the extracted strings and the comments:
            # 1) split: the string is split when the specified character is encountered; 
            # 2) we select the portion of string with index inside square brackets;
            # 3) strip: removes leading or trailing whitespaces from the string. 
            sh_vel_x   =   sh_vel_x.split('=')[-1].strip()
            sh_vel_tot = sh_vel_tot.split('=')[-1].strip()
            
            # Convert to needed variable type (float)
            sh_vel_x   = np.float64(sh_vel_x)
            sh_vel_tot = np.float64(sh_vel_tot)
            
            
            # Read scalar mean stats averaged with different flow realizations.
            # These are time-dependent statistics.
            if numscalar == 1:

                M = np.loadtxt(f'data_post_te/scalar/mean_stats_scalar_realiz-ts{ts}.txt', skiprows=5, delimiter=',', dtype=np.float64)
                        
                mean_phi  = M[:,0]
                var_phi   = M[:,1]
                mean_uphi = M[:,2]
                mean_vphi = M[:,3]
                mean_wphi = M[:,4]
            
            # Opening of 'mean_stats_scalar_realiz-ts' file
            with open(f'data_post_te/scalar/mean_stats_scalar_realiz-ts{ts}.txt', 'r') as file:
    
                # Read all lines into a list
                lines = file.readlines()
            
                # Mean scalar gradient at the wall reading: as always, index is 1 less of the line number (Python convention)
                mg_phi_w = lines[2]   # streamwise shear velocity (based on streamwise mean gradient)  

                # Removing characters in front of the extracted strings and the comments:
                # 1) split: the string is split when the specified character is encountered; 
                # 2) we select the portion of string with index inside square brackets;
                # 3) strip: removes leading or trailing whitespaces from the string. 
                mg_phi_w = mg_phi_w.split('=')[-1].strip()
            
                # Convert to needed variable type (float)
                mg_phi_w = np.float64(mg_phi_w)
            
                # Calculate phi_tau (counterpart of u_tau for the scalar field) (Kozul et al. (2016))
                # We are assuming unitary molecular Prandtl number (thus we are multiplying by kinematic viscosity) 
                # To be verified the correct rescaling for wall oscillations (u_taux or u_tau)
                phi_tau = nu * mg_phi_w / sh_vel_tot
    
    
    # Shear velocities calculation for statistics obtained from 'post_incompact3d', valid for both Channel and TTBL
    if i_switch_plot == False:
    
        """           
        (Total) wall shear stress is used to check maximum numerical resolutions (mesh spacings and viscous time).
        In case of fixed wall(s), the mean spanwise gradient 'mg_z' is zero, so the total wall shear stress 
        is equivalent to the streamwise wall shear stress.
        """
        tau_wtot  = nu*np.sqrt(mg_x[0]**2 + mg_z[0]**2) 

        """
        Streamwise wall shear stress is used to rescale statistics in wall units 
        (with or without wall oscillations).
        """
        tau_wx   = nu*np.abs(mg_x[0])

        # Shear velocities
        sh_vel_x   = np.sqrt(tau_wx)     # streamwise shear velocity (based on streamwise mean gradient)  
        sh_vel_tot = np.sqrt(tau_wtot)   # total shear velocity (based on total mean gradient)
        
        # One scalar field
        if numscalar == 1:    
        
            # Scalar friction quantity phi_tau (counterpart of u_tau) (see Kozul et al. (2016) for example)
            # To be verified the correct rescaling for wall oscillations (u_taux or u_tau)
            phi_tau = nu * (mg_phi[0] / sh_vel_tot) 
                 
    print()
 
    return (
    mean_u, mean_w, var_u, var_v, var_w, mean_uv, 
    vort_x, vort_y, vort_z, mg_x, mg_z, mg_phi, 
    eps, Ruuz, Rvvz, Rwwz, Ruvz, Rssz,
    tke_turbt, tke_presst, tke_difft, tke_prod, tke_pseps,
    snap_numb, i_switch_plot, ts, sh_vel_x, sh_vel_tot, phi_tau
    )

#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Python function to read reference data.  
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def read_ref_data():

    #!--- Initialize reference data ---!
    
    # Lee & Moser (2015), Channel Re_tau = 180
    y_plus_lm              = 0.0
    mean_u_lm              = 0.0
    var_u_lm               = 0.0
    var_v_lm               = 0.0
    var_w_lm               = 0.0
    mean_uv_lm             = 0.0
    
    # Kim et al. (1987), Channel Re_tau = 180
    rz_plus_cuuz_kim       = 0.0
    cuuz_kim               = 0.0
    rz_plus_cvvz_kim       = 0.0
    cvvz_kim               = 0.0
    rz_plus_cwwz_kim       = 0.0
    cwwz_kim               = 0.0
    
    # Touber & Leschziner (2012)
    y_plus_touber          = 0.0
    mean_u_touber          = 0.0
    
    # Yao et al. (2019)
    y_plus_umean_yao       = 0.0
    mean_u_yao             = 0.0
    y_plus_uvar_yao        = 0.0
    var_u_yao              = 0.0
    y_plus_vvar_yao        = 0.0
    var_v_yao              = 0.0
    y_plus_wvar_yao        = 0.0
    var_w_yao              = 0.0
    y_plus_uvmean_yao      = 0.0
    mean_uv_yao            = 0.0
    
    # Moser et al. (1999)
    y_plus_moser_1999      = 0.0
    p_eps_ratio_moser_1999 = 0.0  # ratio of production and dissipation of TKE
    
    # Lee & Moser (2015), Channel Re_tau = 1000
    y_plus_lm1000          = 0.0
    p_eps_ratio_lm1000     = 0.0  # ratio of production and dissipation of TKE minus 1.0
    
    # Kozul et al. (2016)
    y_plus_umean_kozul     = 0.0
    mean_u_kozul           = 0.0
    y_plus_uvar_kozul      = 0.0
    var_u_kozul            = 0.0
    y_plus_vvar_kozul      = 0.0
    var_v_kozul            = 0.0
    y_plus_uvmean_kozul    = 0.0
    mean_uv_kozul          = 0.0
    
    # Mansour et al. (1988)
    y_plus_tke_turbt_mansour  = 0.0  
    tke_turbt_mansour         = 0.0  
    y_plus_tke_presst_mansour = 0.0
    tke_presst_mansour        = 0.0
    y_plus_tke_difft_mansour  = 0.0
    tke_difft_mansour         = 0.0
    y_plus_tke_prod_mansour   = 0.0
    tke_prod_mansour          = 0.0
    y_plus_tke_pseps_mansour  = 0.0
    tke_pseps_mansour         = 0.0
         
    #!--- Reading of fixed walls channels data ---!
                  
    # Reading of Lee & Moser (2015) data
    M = np.loadtxt(dirname + '/lee&moser2015/mean_stats_lee&moser2015.txt', skiprows=72, dtype=np.float64)
    y_plus_lm = M[:,1]
    mean_u_lm = M[:,2]
    
    M = np.loadtxt(dirname + '/lee&moser2015/var_stats_lee&moser2015.txt', skiprows=75, dtype=np.float64)
    var_u_lm   =   M[:,2]
    var_v_lm   =   M[:,3]
    var_w_lm   =   M[:,4]
    mean_uv_lm = - M[:,5]
    
    # Velocity auto-correlations, Kim et al. (1987) data, y+ = 10.52
    M = np.loadtxt(dirname + '/kim1987/cuuz_kim1987.txt', skiprows=7, delimiter=',', dtype=np.float64)
    rz_plus_cuuz_kim = M[:,0]
    cuuz_kim         = M[:,1]
    
    M = np.loadtxt(dirname + '/kim1987/cvvz_kim1987.txt', skiprows=7, delimiter=',', dtype=np.float64)
    rz_plus_cvvz_kim = M[:,0]
    cvvz_kim         = M[:,1]
    
    M = np.loadtxt(dirname + '/kim1987/cwwz_kim1987.txt', skiprows=7, delimiter=',', dtype=np.float64)
    rz_plus_cwwz_kim = M[:,0]
    cwwz_kim         = M[:,1]
        
    #!--- Reading of oscillating walls channels data (A^+ = 12, T^+ = 100, Re_tau = 200) ---! 
      
    # Mean velocity profile (Touber & Leschziner (2012))
    M = np.loadtxt(dirname + '/touber2012/umean_touber2012.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_touber = M[:,0]
    mean_u_touber = M[:,1]
        
    # Mean velocity profile, Reynolds stress and RMSs (Yao et al. (2019))
    M = np.loadtxt(dirname + '/yao2019/umean_yao2019.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_umean_yao = M[:,0]
    mean_u_yao       = M[:,1]
        
    M = np.loadtxt(dirname + '/yao2019/uvar_yao2019.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_uvar_yao = M[:,0]
    var_u_yao       = M[:,1]
        
    M = np.loadtxt(dirname + '/yao2019/vvar_yao2019.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_vvar_yao = M[:,0]
    var_v_yao       = M[:,1]
                
    M = np.loadtxt(dirname + '/yao2019/vvar_yao2019.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_wvar_yao = M[:,0]
    var_w_yao       = M[:,1]
        
    M = np.loadtxt(dirname + '/yao2019/uvmean_yao2019.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_uvmean_yao = M[:,0]
    mean_uv_yao       = M[:,1]
        
    # Rescale Yao et al. (2019) data from uncontrolled to controlled shear velocity (0: uncontrolled, c: controlled)
    cf_0_yao = 0.00790
    cf_c_yao = 0.00511
        
    sh_vel_0_yao = (2.0/3.0)*np.sqrt(cf_0_yao / 2.0)
    sh_vel_c_yao = (2.0/3.0)*np.sqrt(cf_c_yao / 2.0)
        
    # Mean velocity is already rescaled by the actual shear velocity in Yao's data
        
    # Rescale RMSs and obtain variances
    y_plus_uvar_yao   = y_plus_uvar_yao   * (sh_vel_c_yao / sh_vel_0_yao)
    y_plus_vvar_yao   = y_plus_vvar_yao   * (sh_vel_c_yao / sh_vel_0_yao)
    y_plus_wvar_yao   = y_plus_wvar_yao   * (sh_vel_c_yao / sh_vel_0_yao)
    y_plus_uvmean_yao = y_plus_uvmean_yao * (sh_vel_c_yao / sh_vel_0_yao)
        
    var_u_yao   = (var_u_yao   *  sh_vel_0_yao / sh_vel_c_yao)**2
    var_v_yao   = (var_v_yao   *  sh_vel_0_yao / sh_vel_c_yao)**2
    var_w_yao   = (var_w_yao   *  sh_vel_0_yao / sh_vel_c_yao)**2
    mean_uv_yao = mean_uv_yao  * (sh_vel_0_yao / sh_vel_c_yao)**2
    
    # Moser et al. (1999)
    M = np.loadtxt(dirname + '/moser1999/p_eps_ratio_moser1999.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_moser_1999      = M[:,0]
    p_eps_ratio_moser_1999 = M[:,1]  
    
    # Lee & Moser (2015), Channel Re_tau = 1000
    M = np.loadtxt(dirname + '/lee&moser2015/p_eps_ratio_minus1_lee&moser2015.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_lm1000      = M[:,0]
    p_eps_ratio_lm1000 = M[:,1]
    
    #!--- Kozul et al. (2016) ---!
    M = np.loadtxt(dirname + '/kozul2016/umean_kozul2016.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_umean_kozul     = M[:,0]
    mean_u_kozul           = M[:,1]
    
    M = np.loadtxt(dirname + '/kozul2016/uprimerms_kozul2016.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_uvar_kozul      = M[:,0]
    var_u_kozul            = M[:,1]
    
    M = np.loadtxt(dirname + '/kozul2016/vprimerms_kozul2016.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_vvar_kozul      = M[:,0]
    var_v_kozul            = M[:,1]
    
    M = np.loadtxt(dirname + '/kozul2016/uvmean_kozul2016.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_uvmean_kozul    = M[:,0]
    mean_uv_kozul          = M[:,1]
     
    # Square RMSs to obtain variances
    var_u_kozul = var_u_kozul**2
    var_v_kozul = var_v_kozul**2
    
    #!--- Mansour et al. (1988) ---!
    M = np.loadtxt(dirname + '/mansour1988/tke_turbt_mansour1988.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_tke_turbt_mansour  = M[:,0]
    tke_turbt_mansour         = M[:,1]
    
    M = np.loadtxt(dirname + '/mansour1988/tke_presst_mansour1988.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_tke_presst_mansour = M[:,0]
    tke_presst_mansour        = M[:,1]
    
    M = np.loadtxt(dirname + '/mansour1988/tke_difft_mansour1988.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_tke_difft_mansour  = M[:,0]
    tke_difft_mansour         = M[:,1]
    
    M = np.loadtxt(dirname + '/mansour1988/tke_prod_mansour1988.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_tke_prod_mansour   = M[:,0]
    tke_prod_mansour          = M[:,1]
    
    M = np.loadtxt(dirname + '/mansour1988/tke_pseps_mansour1988.txt', skiprows=8, delimiter=',', dtype=np.float64)
    y_plus_tke_pseps_mansour  = M[:,0]
    tke_pseps_mansour         = M[:,1]
          
    # Return to main program with extracted reference data
    return (
    y_plus_lm,                 mean_u_lm, var_u_lm, var_v_lm, var_w_lm, mean_uv_lm,
    rz_plus_cuuz_kim,          cuuz_kim, 
    rz_plus_cvvz_kim,          cvvz_kim,
    rz_plus_cwwz_kim,          cwwz_kim,
    y_plus_touber,             mean_u_touber,
    y_plus_umean_yao,          mean_u_yao,
    y_plus_uvar_yao,           var_u_yao,
    y_plus_vvar_yao,           var_v_yao,
    y_plus_wvar_yao,           var_w_yao,
    y_plus_uvmean_yao,         mean_uv_yao, 
    y_plus_moser_1999,         p_eps_ratio_moser_1999,
    y_plus_lm1000,             p_eps_ratio_lm1000,
    y_plus_umean_kozul,        mean_u_kozul,
    y_plus_uvar_kozul,         var_u_kozul,
    y_plus_vvar_kozul,         var_v_kozul,
    y_plus_uvmean_kozul,       mean_uv_kozul,
    y_plus_tke_turbt_mansour,  tke_turbt_mansour,  
    y_plus_tke_presst_mansour, tke_presst_mansour,
    y_plus_tke_difft_mansour,  tke_difft_mansour,
    y_plus_tke_prod_mansour,   tke_prod_mansour,
    y_plus_tke_pseps_mansour,  tke_pseps_mansour    
    ) 

#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!      
  
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Python function to read reference evolution quantities 
!              for TTBL. Reference data from G. Boga
!              (Cimarelli et al. (2024a)).  
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def read_ref_data_temp_evol():

    #!--- Initialize reference data ---!
    
    # G. Boga
    t_gboga             = 0.0
    retau_vs_time_gboga = 0.0
    
    retau_gboga         = 0.0
    retheta_gboga       = 0.0
    utau_gboga          = 0.0
    delta_99_gboga      = 0.0
    disp_t_gboga        = 0.0
    cf_gboga            = 0.0
        
    # Reading of G. Boga data (Cimarelli et al. (2024a))
    M = np.loadtxt(dirname + '/gboga/retau500/Re_tau500_time_largettblnr4_boga.dat', skiprows=1, dtype=np.float64)
    t_gboga             = M[:,0]
    retau_vs_time_gboga = M[:,1]
    
    M = np.loadtxt(dirname + '/gboga/retau500/temp_evolution500_largettblnr4_boga.dat', skiprows=1, dtype=np.float64)
    retau_gboga    = M[:,0]
    retheta_gboga  = M[:,1]  
    utau_gboga     = M[:,2]
    delta_99_gboga = M[:,3]
    disp_t_gboga   = M[:,4]
    cf_gboga       = M[:,6]
        
    # Return to main program with extracted reference data
    return (
    t_gboga, 
    retau_vs_time_gboga,
    retau_gboga,
    retheta_gboga,
    utau_gboga,
    delta_99_gboga,
    disp_t_gboga,
    cf_gboga    
    ) 

#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!      
      
    
