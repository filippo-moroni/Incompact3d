#!-----------------------------------------!
#! Python function to read reference data. !
#!-----------------------------------------!

import numpy as np

def read_reference_data(itype, iswitch_wo):

    # Reference data for Channel flow simulations
    if itype == 3:
               
        # Reading of Lee & Moser (2015) data
        M = np.loadtxt('reference_data/lee&moser2015/mean_stats_lee&moser2015.txt', skiprows=72, dtype=np.float64)
        y_plus_lm = M[:,1]
        mean_u_lm = M[:,2]
    
        M = np.loadtxt('reference_data/lee&moser2015/var_stats_lee&moser2015.txt', skiprows=75, dtype=np.float64)
        var_u_lm   =   M[:,2]
        var_v_lm   =   M[:,3]
        mean_uv_lm = - M[:,5]
    
        # Velocity auto-correlations, Kim et al. (1987) data, y+ = 10.52
        M = np.loadtxt('reference_data/kim1987/cuuz_kim1987.txt', skiprows=7, delimiter=',', dtype=np.float64)
        rz_plus_cuuz_kim = M[:,0]
        cuuz_kim         = M[:,1]
    
        M = np.loadtxt('reference_data/kim1987/cvvz_kim1987.txt', skiprows=7, delimiter=',', dtype=np.float64)
        rz_plus_cvvz_kim = M[:,0]
        cvvz_kim         = M[:,1]
    
        M = np.loadtxt('reference_data/kim1987/cwwz_kim1987.txt', skiprows=7, delimiter=',', dtype=np.float64)
        rz_plus_cwwz_kim = M[:,0]
        cwwz_kim         = M[:,1]
        
        # Reading of wall-oscillations data (A^+ = 12, T^+ = 100) 
        if iswitch_wo == 1:
    
            # Mean velocity profile (Touber & Leschziner (2012))
            M = np.loadtxt('reference_data/touber2012/umean_touber2012.txt', skiprows=8, delimiter=',', dtype=np.float64)
            y_plus_touber = M[:,0]
            mean_u_touber = M[:,1]
        
            # Mean velocity profile, Reynolds stress and RMSs (Yao et al. (2019))
            M = np.loadtxt('reference_data/yao2019/umean_yao2019.txt', skiprows=8, delimiter=',', dtype=np.float64)
            y_plus_umean_yao = M[:,0]
            mean_u_yao       = M[:,1]
        
            M = np.loadtxt('reference_data/yao2019/uvar_yao2019.txt', skiprows=8, delimiter=',', dtype=np.float64)
            y_plus_uvar_yao = M[:,0]
            var_u_yao       = M[:,1]
        
            M = np.loadtxt('reference_data/yao2019/vvar_yao2019.txt', skiprows=8, delimiter=',', dtype=np.float64)
            y_plus_vvar_yao = M[:,0]
            var_v_yao       = M[:,1]
        
            M = np.loadtxt('reference_data/yao2019/uvmean_yao2019.txt', skiprows=8, delimiter=',', dtype=np.float64)
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
            y_plus_uvmean_yao = y_plus_uvmean_yao * (sh_vel_c_yao / sh_vel_0_yao)
        
            var_u_yao   = (var_u_yao   *  sh_vel_0_yao / sh_vel_c_yao)**2
            var_v_yao   = (var_v_yao   *  sh_vel_0_yao / sh_vel_c_yao)**2
            mean_uv_yao = mean_uv_yao  * (sh_vel_0_yao / sh_vel_c_yao)**2
    
    # Return to main program with extracted reference data
    return 
    y_plus_lm,         mean_u_lm, var_u_lm, var_v_lm, mean_uv_lm,
    rz_plus_cuuz_kim,  cuuz_kim, 
    rz_plus_cvvz_kim,  cvvz_kim,
    rz_plus_cwwz_kim,  cwwz_kim,
    y_plus_touber,     mean_u_touber,
    y_plus_umean_yao,  mean_u_yao,
    y_plus_uvar_yao,   var_u_yao,
    y_plus_vvar_yao,   var_v_yao,
    y_plus_uvmean_yao, mean_uv_yao         
        
        
