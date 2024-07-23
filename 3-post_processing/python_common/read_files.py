#!------------------------------------------------------------!
#! In this file, we store Python functions to read .txt files !
#! for post-processing and plots:                             !
#! - read_input_files: to read 'input.i3d' and 'post.prm';    !
#! - read_data:        to read statistics data;               !
#! - read_ref_data:    to read reference data.                !
#!------------------------------------------------------------!

#!-------------------------------------!
#! Function to read Incompact3d files: !
#! 'input.i3d' and 'post.prm'.         !
#!-------------------------------------!

import numpy as np

# Define 'read_input_files' function, that reads parameters of simulation 
# from 'input.i3d' and 'post.prm' files.

def read_input_files(filename1,filename2):
      
    # Opening of 'input.i3d' file
    with open(filename1, 'r') as file:
        
        # Read all lines into a list
        lines = file.readlines()
    
        # Extract itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, Re, dt, numscalar, iswitch_wo 
        # As always, index is 1 less of the line number (Python convention)
        itype      = lines[7]  
        nx         = lines[14]
        ny         = lines[15]
        nz         = lines[16]
        istret     = lines[17]
        beta       = lines[18]
        Lx         = lines[21]
        Ly         = lines[22]
        Lz         = lines[23]
        re         = lines[26]
        dt         = lines[29]
        numscalar  = lines[34]
        iswitch_wo = lines[88]
    
        # Removing characters in front of the extracted strings and the comments:
        # 1) split: the string is split when the specified character is encountered; 
        # 2) we select the portion of string with index inside square brackets;
        # 3) strip: removes leading or trailing whitespaces from the string. 
    
        itype = itype.split('=')[-1].strip()
    
        nx         = nx.split('!')[0]
        nx         = nx.split('=')[-1].strip()
        
        ny         = ny.split('!')[0]
        ny         = ny.split('=')[-1].strip()
    
        nz         = nz.split('!')[0]
        nz         = nz.split('=')[-1].strip()
        
        istret     = istret.split('!')[0]
        istret     = istret.split('=')[-1].strip()
        
        beta       = beta.split('!')[0]
        beta       = beta.split('=')[-1].strip()
    
        Lx         = Lx.split('!')[0]
        Lx         = Lx.split('=')[-1].strip()
    
        Ly         = Ly.split('!')[0]
        Ly         = Ly.split('=')[-1].strip()
    
        Lz         = Lz.split('!')[0]
        Lz         = Lz.split('=')[-1].strip()
    
        re         = re.split('!')[0]
        re         = re.split('=')[-1].strip()
        
        dt         = dt.split('!')[0]
        dt         = dt.split('=')[-1].strip()
        
        numscalar  = numscalar.split('!')[0]
        numscalar  = numscalar.split('=')[-1].strip()
    
        iswitch_wo = iswitch_wo.split('!')[0]
        iswitch_wo = iswitch_wo.split('=')[-1].strip()
    
        # Convert to needed variable type (integer, float, etc.)
        itype      = int(itype)
        nx         = int(nx)
        ny         = int(ny)
        nz         = int(nz)
        istret     = int(istret)
        beta       = np.float64(beta)
        Lx         = np.float64(Lx)
        Ly         = np.float64(Ly)
        Lz         = np.float64(Lz)
        re         = np.float64(re)
        dt         = np.float64(dt)
        numscalar  = int(numscalar)
        iswitch_wo = int(iswitch_wo)
    
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
    return itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, numscalar, iswitch_wo, add_string, file1, filen, icrfile, nr, post_mean, post_vort, post_diss, post_corz, post_tke_eq 

#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
    
#!------------------------------------------!
#! Python function to read statistics data. !
#!------------------------------------------!

import numpy as np

def read_data(itype, numscalar, post_mean, post_vort, post_diss, post_corz, post_tke_eq, snap_numb=None):

    # Initialize variables
    snap_numb   = None
    mean_u      = 0.0
    mean_w      = 0.0
    var_u       = 0.0
    var_v       = 0.0
    var_w       = 0.0
    mean_uv     = 0.0
    
    vort_x      = 0.0
    vort_y      = 0.0
    vort_z      = 0.0
    mg_tot      = 0.0
    mg_x        = 0.0
    mg_z        = 0.0
    
    eps         = 0.0
    
    Ruuz        = 0.0
    Rvvz        = 0.0
    Rwwz        = 0.0
    Ruvz        = 0.0
    Rppz        = 0.0
    
    tke_conv    = 0.0
    tke_turbt   = 0.0   
    tke_pstrain = 0.0     
    tke_difft   = 0.0
    tke_prod    = 0.0
    tke_pseps   = 0.0
    
    print()

    # Channel
    if itype == 3:
    
        print("!--- Plotting of statistics for a channel ---!")

        # Setting 'snap_numb' to an empty string
        snap_numb = ''

        # Reading of mean statistics
        if post_mean:
            M = np.loadtxt('data_post/mean_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
            mean_u  = M[:,0]
            mean_w  = M[:,2]   
            var_u   = M[:,3]
            var_v   = M[:,4]
            var_w   = M[:,5]
            mean_uv = M[:,12]
    
        # Reading of vorticity components and mean gradient
        if post_vort:
            M = np.loadtxt('data_post/vort_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
            vort_x = M[:,0]
            vort_y = M[:,1]
            vort_z = M[:,2]
            mg_tot = M[:,3]
            mg_x   = M[:,4]
            mg_z   = M[:,5]
    
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
                
                Rppz = np.loadtxt('data_post/Rppz.txt', skiprows=0, delimiter=None, dtype=np.float64)
        
        # Read of TKE 
        if post_tke_eq:
            M = np.loadtxt('data_post/tke_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
            tke_conv    = M[:,0]
            tke_turbt   = M[:,1]   
            tke_pstrain = M[:,2]     
            tke_difft   = M[:,3]
            tke_prod    = M[:,4]
            tke_pseps   = M[:,5]
        
    # TTBL
    elif itype == 13:

        print("!--- Plotting of statistics for a TTBL ---!")
        print()
    
        # Asking to the user the specific snapshot to show
        snap_numb = input("Enter the snapshot number to show: ")
    
        # Pad with zeros to match snapshots' naming  
        snap_numb = snap_numb.zfill(4)
             
        # Reading of mean statistics
        if post_mean:
            M = np.loadtxt(f'data_post/mean_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
            mean_u  = M[:,0]
            mean_w  = M[:,2]   
            var_u   = M[:,3]
            var_v   = M[:,4]
            var_w   = M[:,5]
            mean_uv = M[:,12]
    
        # Reading of vorticity components and mean gradient
        if post_vort:
            M = np.loadtxt(f'data_post/vort_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
            vort_x = M[:,0]
            vort_y = M[:,1]
            vort_z = M[:,2]
            mg_tot = M[:,3]
            mg_x   = M[:,4]
            mg_z   = M[:,5]
    
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
            
                Rppz = np.loadtxt(f'data_post/Rppz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)
        
        # Read of TKE
        if post_tke_eq: 
            M = np.loadtxt(f'data_post/tke_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
            tke_conv    = M[:,0]
            tke_turbt   = M[:,1]   
            tke_pstrain = M[:,2]     
            tke_difft   = M[:,3]
            tke_prod    = M[:,4]
            tke_pseps   = M[:,5]

    print()
 
    return (
    mean_u, mean_w, var_u, var_v, var_w, mean_uv, 
    vort_x, vort_y, vort_z, mg_tot, mg_x, mg_z,
    eps, Ruuz, Rvvz, Rwwz, Ruvz, Rppz,
    tke_conv, tke_turbt, tke_pstrain, tke_difft, tke_prod, tke_pseps,
    snap_numb
    )

#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    

#!-----------------------------------------!
#! Python function to read reference data. !
#!-----------------------------------------!

import numpy as np

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
    
    #!--- Reading of fixed walls channels data ---!
    
    # Directory where reference data are stored
    dirname = '/home/n286654/Desktop/Incompact3d/3-post_processing/reference_data'
              
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
    
    # Return to main program with extracted reference data
    return (
    y_plus_lm,          mean_u_lm, 
    var_u_lm, var_v_lm, var_w_lm,      mean_uv_lm,
    rz_plus_cuuz_kim,   cuuz_kim, 
    rz_plus_cvvz_kim,   cvvz_kim,
    rz_plus_cwwz_kim,   cwwz_kim,
    y_plus_touber,      mean_u_touber,
    y_plus_umean_yao,   mean_u_yao,
    y_plus_uvar_yao,    var_u_yao,
    y_plus_vvar_yao,    var_v_yao,
    y_plus_wvar_yao,    var_w_yao,
    y_plus_uvmean_yao,  mean_uv_yao, 
    y_plus_moser_1999,  p_eps_ratio_moser_1999,
    y_plus_lm1000,      p_eps_ratio_lm1000
    ) 

#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!      
  
#!------------------------------------------------------------------!
#! Python function to read reference evolution quantities for TTBL. !
#!------------------------------------------------------------------!

import numpy as np

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
    
    # Directory where reference data are stored
    dirname = '/home/n286654/Desktop/Incompact3d/3-post_processing/reference_data'
    
    # Reading of G. Boga data
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
      
    
