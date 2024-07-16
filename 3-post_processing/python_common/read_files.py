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
    
        # Extract itype, nx, ny, nz, Lx, Ly, Lz, Re, numscalar, iswitch_wo 
        # As always, index is 1 less of the line number (Python convention)
        itype      = lines[7]  
        nx         = lines[14]
        ny         = lines[15]
        nz         = lines[16]
        Lx         = lines[21]
        Ly         = lines[22]
        Lz         = lines[23]
        re         = lines[26]
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
    
        Lx         = Lx.split('!')[0]
        Lx         = Lx.split('=')[-1].strip()
    
        Ly         = Ly.split('!')[0]
        Ly         = Ly.split('=')[-1].strip()
    
        Lz         = Lz.split('!')[0]
        Lz         = Lz.split('=')[-1].strip()
    
        re         = re.split('!')[0]
        re         = re.split('=')[-1].strip()
        
        numscalar  = numscalar.split('!')[0]
        numscalar  = numscalar.split('=')[-1].strip()
    
        iswitch_wo = iswitch_wo.split('!')[0]
        iswitch_wo = iswitch_wo.split('=')[-1].strip()
    
        # Convert to needed variable type (integer, float, etc.)
        itype      = int(itype)
        nx         = int(nx)
        ny         = int(ny)
        nz         = int(nz)
        Lx         = np.float64(Lx)
        Ly         = np.float64(Ly)
        Lz         = np.float64(Lz)
        re         = np.float64(re)
        numscalar  = int(numscalar)
        iswitch_wo = int(iswitch_wo)
    
    # Opening of 'post.prm' file
    with open(filename2, 'r') as file:
        
        # Read all lines into a list
        lines = file.readlines()
   
        # Extract needed lines  
        add_string = lines[3]   # Flow case name
        file1      = lines[5]   # First snapshot index
        filen      = lines[6]   # Final snapshot index
        icrfile    = lines[7]   # File increment
        nr         = lines[8]   # Number of flow realizations
        
        # Extract the needed variables
        file1   =   file1.split('#')[0].strip()
        filen   =   filen.split('#')[0].strip()
        icrfile = icrfile.split('#')[0].strip()
        nr      =      nr.split('#')[0].strip()
        
        add_string = add_string.split('!')[0]
        add_string = add_string.rstrip()
        
        # Convert to needed variable type (integer)
        file1      = int(file1)
        filen      = int(filen)
        icrfile    = int(icrfile)
        nr         = int(nr)
    
    # Halve the points in y direction for a channel
    if itype == 3:
        ny = (ny - 1) // 2 + 1
        
    # Return to main program with extracted parameters
    return itype, nx, ny, nz, Lx, Ly, Lz, re, numscalar, iswitch_wo, file1, filen, icrfile, nr, add_string

#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
    
#!------------------------------------------!
#! Python function to read statistics data. !
#!------------------------------------------!

import numpy as np

def read_data(itype, numscalar):

    # Initialize variables
    snap_numb = ''
    mean_u    = 0.0
    mean_w    = 0.0
    var_u     = 0.0
    var_v     = 0.0
    mean_uv   = 0.0
    
    vort_x    = 0.0
    vort_y    = 0.0
    vort_z    = 0.0
    mg_tot    = 0.0
    mg_x      = 0.0
    mg_z      = 0.0
    
    Ruuz      = 0.0
    Rvvz      = 0.0
    Rwwz      = 0.0
    Ruvz      = 0.0
    Rppz      = 0.0
    
    print()

    # Channel
    if itype == 3:
    
        print("!--- Plotting of statistics for a channel ---!")

        # Reading of mean statistics
        M1 = np.loadtxt('data_post/mean_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
        # Reading of vorticity components and mean gradient
        M2 = np.loadtxt('data_post/vort_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
        # Reading of the mean total dissipation
        eps = np.loadtxt('data_post/diss_stats.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
        # Reading of correlations
        Ruuz = np.loadtxt('data_post/Ruuz.txt', skiprows=0, delimiter=None, dtype=np.float64)
        Rvvz = np.loadtxt('data_post/Rvvz.txt', skiprows=0, delimiter=None, dtype=np.float64)
        Rwwz = np.loadtxt('data_post/Rwwz.txt', skiprows=0, delimiter=None, dtype=np.float64)
        Ruvz = np.loadtxt('data_post/Ruvz.txt', skiprows=0, delimiter=None, dtype=np.float64)
      
        # Read scalar field statistics
        if numscalar == 1:
                
            Rppz = np.loadtxt('data_post/Rppz.txt', skiprows=0, delimiter=None, dtype=np.float64)
        
    # TTBL
    elif itype == 13:

        print("!--- Plotting of statistics for a TTBL ---!")
        print()
    
        # Asking to the user the specific snapshot to show
        snap_numb = input("Enter the snapshot number to show: ")
    
        # Pad with zeros to match snapshots' naming  
        snap_numb = snap_numb.zfill(4)
             
        # Reading of mean statistics
        M1 = np.loadtxt(f'data_post/mean_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
        # Reading of vorticity components and mean gradient
        M2 = np.loadtxt(f'data_post/vort_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
    
        # Reading of the mean total dissipation
        eps = np.loadtxt(f'data_post/diss_stats-{snap_numb}.txt', skiprows=1, delimiter=',', dtype=np.float64)
        
        # Reading of correlations
        Ruuz = np.loadtxt(f'data_post/Ruuz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)
        Rvvz = np.loadtxt(f'data_post/Rvvz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)
        Rwwz = np.loadtxt(f'data_post/Rwwz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)
        Ruvz = np.loadtxt(f'data_post/Ruvz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)
        
        # Read scalar field statistics
        if numscalar == 1:
            
            Rppz = np.loadtxt(f'data_post/Rppz-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)

    print()

    # Extracting quantities from the full matrices
    mean_u  = M1[:,0]
    mean_w  = M1[:,2]   
    var_u   = M1[:,3]
    var_v   = M1[:,4]
    mean_uv = M1[:,12]

    vort_x = M2[:,0]
    vort_y = M2[:,1]
    vort_z = M2[:,2]
    mg_tot = M2[:,3]
    mg_x   = M2[:,4]
    mg_z   = M2[:,5]
    
    return (
    mean_u, mean_w, var_u, var_v, mean_uv, 
    vort_x, vort_y, vort_z, mg_tot, mg_x, mg_z,
    eps, Ruuz, Rvvz, Rwwz, Ruvz, Rppz,
    snap_numb
    )

#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    

#!-----------------------------------------!
#! Python function to read reference data. !
#!-----------------------------------------!

import numpy as np

def read_ref_data():

    # Initialize reference data
    y_plus_lm         = 0.0
    mean_u_lm         = 0.0
    var_u_lm          = 0.0
    var_v_lm          = 0.0
    mean_uv_lm        = 0.0
    rz_plus_cuuz_kim  = 0.0
    cuuz_kim          = 0.0
    rz_plus_cvvz_kim  = 0.0
    cvvz_kim          = 0.0
    rz_plus_cwwz_kim  = 0.0
    cwwz_kim          = 0.0
    y_plus_touber     = 0.0
    mean_u_touber     = 0.0
    y_plus_umean_yao  = 0.0
    mean_u_yao        = 0.0
    y_plus_uvar_yao   = 0.0
    var_u_yao         = 0.0
    y_plus_vvar_yao   = 0.0
    var_v_yao         = 0.0
    y_plus_uvmean_yao = 0.0
    mean_uv_yao       = 0.0  
    
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
        
    #!--- Reading of oscillating walls channels data (A^+ = 12, T^+ = 100) ---! 
      
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
    y_plus_uvmean_yao = y_plus_uvmean_yao * (sh_vel_c_yao / sh_vel_0_yao)
        
    var_u_yao   = (var_u_yao   *  sh_vel_0_yao / sh_vel_c_yao)**2
    var_v_yao   = (var_v_yao   *  sh_vel_0_yao / sh_vel_c_yao)**2
    mean_uv_yao = mean_uv_yao  * (sh_vel_0_yao / sh_vel_c_yao)**2
    
    # Return to main program with extracted reference data
    return (
    y_plus_lm,         mean_u_lm, var_u_lm, var_v_lm, mean_uv_lm,
    rz_plus_cuuz_kim,  cuuz_kim, 
    rz_plus_cvvz_kim,  cvvz_kim,
    rz_plus_cwwz_kim,  cwwz_kim,
    y_plus_touber,     mean_u_touber,
    y_plus_umean_yao,  mean_u_yao,
    y_plus_uvar_yao,   var_u_yao,
    y_plus_vvar_yao,   var_v_yao,
    y_plus_uvmean_yao, mean_uv_yao 
    ) 

#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!      
  
    
    
