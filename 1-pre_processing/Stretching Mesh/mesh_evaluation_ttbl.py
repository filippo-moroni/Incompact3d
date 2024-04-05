#!---------------------------------------------------------!
#! mesh_evaluation_ttbl.py:                                !
#! improved Python version of 'mesh_evaluation.f90'.       !
#!                                                         !
#! Info:                                                   !
#!  - Same values are obtained from both codes (tested).   !
#!  - y-coordinates at the center of mesh nodes are not    !
#!    evaluated here.                                      !
#!---------------------------------------------------------!

#!------------------- Description: ------------------------!
# 
# In this script we calculate mesh, numerics and 
# flow parameters for temporal TBL simulations
# (see the related .pdf file for more details).
#!---------------------------------------------------------!

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
from tabulate import tabulate

# Font settings
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

# Inputs
istret = 3               # y mesh refinement (0:no, 1:center, 2:both sides, 3:bottom)
beta = 3.0               # beta parameter for mesh stretching
nu = 0.002               # kinematic viscosity (if D = 1 and U_wall = 1, Re_D = 500)
uwall = 1.0              # velocity of the wall
delta_t = 0.01           # time-step
twd = 1.0                # trip wire diameter D
re = 1.0/nu              # Reynolds number as defined in Incompact3d

bl_thickness = 21.7*twd  # displacement thickness of the temporal TBL at Re_tau = 500 (Cimarelli et al. (2024))
cf = 0.007               # maximum cf estimated at peak (Cimarelli et al. (2024))

nx = 34                  # number of points in x direction
ny = 651                 # number of points in y direction
nz = 52                  # number of points in z direction

xlx = 10.0               # domain dimension in x direction
yly = 3.0*bl_thickness   # domain dimension in y direction
zlz = 8.0                # domain dimension in z direction

nym = ny - 1             # if periodic BC is imposed along y, nym = ny, otherwise nym = ny - 1

# Declare local constant Pi
pi = np.pi

# Work variables
yeta = np.zeros(ny)
yp = np.zeros(ny)
Uo = np.zeros(ny)

# Total number of points
ntot = nx*ny*nz

# Start of calculations
yinf = -yly/2.0
den = 2.0 * beta * yinf
xnum = -yinf - np.sqrt(pi*pi*beta*beta + yinf*yinf)
alpha = np.abs(xnum/den)
xcx = 1.0/beta/alpha

if alpha != 0.0:
    for j in range(1, ny):
        if istret == 1:
            yeta[j] = j*(1.0/nym)
        elif istret == 2:
            yeta[j] = j*(1.0/nym) - 0.5
        elif istret == 3:
            yeta[j] = j*(0.5/nym) - 0.5

        den1 = np.sqrt(alpha*beta + 1.0)
        xnum = den1/np.sqrt(alpha/pi)/np.sqrt(beta)/np.sqrt(pi)
        den = 2.0 * np.sqrt(alpha/pi) * np.sqrt(beta) * pi * np.sqrt(pi)
        den3 = ((np.sin(pi*yeta[j]))**2/beta/pi) + alpha/pi
        den4 = 2.0*alpha*beta - np.cos(2.0*pi*yeta[j]) + 1.0
        xnum1 = (np.arctan(xnum*np.tan(pi*yeta[j]))) * den4 / den1 / den3 / den
        cst = np.sqrt(beta)*pi / (2.0*np.sqrt(alpha)*np.sqrt(alpha*beta + 1.0))

        if istret == 1:
            if yeta[j] < 0.5:
                yp[j] = xnum1 - cst - yinf
            elif yeta[j] == 0.5:
                yp[j] = 0.0 - yinf
            elif yeta[j] > 0.5:
                yp[j] = xnum1 + cst - yinf

        if istret == 2:
            if yeta[j] < 0.5:
                yp[j] = xnum1 - cst + yly
            elif yeta[j] == 0.5:
                yp[j] = 0.0 + yly
            elif yeta[j] > 0.5:
                yp[j] = xnum1 + cst + yly

        if istret == 3:
            if yeta[j] < 0.5:
                yp[j] = (xnum1 - cst + yly) * 2.0
            elif yeta[j] == 0.5:
                yp[j] = (0.0 + yly) * 2.0
            elif yeta[j] > 0.5:
                yp[j] = (xnum1 + cst + yly) * 2.0

if alpha == 0.0:
    yp[0] = -1.e10
    for j in range(1, ny):
        yeta[j] = j*(1.0/ny)
        yp[j] = -beta * np.cos(pi*yeta[j]) / np.sin(yeta[j]*pi)


# Calculate the spacings along x and z (uniform)
delta_x = xlx / nx
delta_z = zlz / nz

# This part is valid for meshes with refinement at the bottom boundary and for temporal TBL cases 
if istret == 3:
   
    # Calculating the initial thickness of the shear layer (see Kozul et al. (2016))
    theta_sl = 54.0 * nu / uwall
    
    # Mean gradient due to initial condition (analytical derivative)  
    mg = - uwall / (4.0 * theta_sl) * (1.0 / np.cosh(twd / 2.0 / theta_sl))**2
    
    
    #!--- Shear velocities (IC, peak cf and Re_tau = 500): ---!
    
    # Shear velocity due to initial condition
    sh_vel_ic = np.sqrt(nu * np.abs(mg))
       
    # Shear velocity peak (reference value of cf = 0.07, Cimarelli et al. (2024))
    sh_vel_peak = np.sqrt((cf/2.0)) * uwall
    
    # Shear velocity at Re_tau = 500:
    sh_vel_500 = 500.0 * nu / bl_thickness
    
    #!--------------------------------------------------------!

    # Calculate viscous unit at IC, at peak cf and at Re_tau = 500
    delta_nu_ic   = nu / sh_vel_ic
    delta_nu_peak = nu / sh_vel_peak
    delta_nu_500  = nu / sh_vel_500

    # First element y-dimension
    delta_y1 = yp[1] - yp[0]
    
    # Rescaling the mesh spacings with viscous unit at peak cf
    delta_x_nd_peak  = delta_x  / delta_nu_peak
    delta_y1_nd_peak = delta_y1 / delta_nu_peak
    delta_z_nd_peak  = delta_z  / delta_nu_peak
        
    # Non-dimensional domain dimensions (nd: non-dimensional) (at IC and at peak cf)  
    xlx_nd_ic = xlx / delta_nu_ic
    yly_nd_ic = yly / delta_nu_ic
    zlz_nd_ic = zlz / delta_nu_ic
    
    xlx_nd_peak = xlx / delta_nu_peak
    yly_nd_peak = yly / delta_nu_peak
    zlz_nd_peak = zlz / delta_nu_peak
    
    #!--- Estimation of numerics-related parameters (CFL, D, Pé, S) at IC ---!
   
    CFL = uwall * delta_t / delta_x      
    D =   nu * delta_t / (delta_x**2)    
    Pe =  uwall * delta_x / nu
        
    # Stability parameter (S < 1) (see Thompson et al. (1985)) 
    S = ((uwall**2)*delta_t)/(2.0*nu)
    
    #!-----------------------------------------------------------------------!
    
    #!--- Check on the initial velocity profile ---!
        
    # Initial velocity profile (tanh) (Kozul et al. (2016))
    for j in range(0, ny):
    	Uo[j] = uwall * (0.5 + 0.5 * (math.tanh((twd/2.0/theta_sl)*(1.0 - yp[j]/twd))))
    
    # Rescaling the initial velocity profile and the y-coordinates
    Uo = Uo / sh_vel_ic
    yp_ic = yp / delta_nu_ic
        
    # Plotting of the initial velocity profile in wall units, first 50 points
    plt.scatter(yp_ic[0:50], Uo[0:50])
    plt.title("Initial velocity profile near the wall", fontsize=30)
    plt.xlabel("$y^+$", fontsize=30)
    plt.ylabel("$U_o^+$", fontsize=30)
    plt.show()
        
    # Calculate the thickness delta99^+ of the initial shear layer
    j = 0
    while j <= ny - 1 and Uo[j] > Uo[0]*0.01:
    	sl_99_ic = yp_ic[j]
    	j = j + 1
    	
    # Calculation of the number of mesh nodes in the initial shear layer
    npsl = 0      # number of points shear layer
    height = 0.0  # cumulative height in viscous unit (y+)
  
    for j in range(1, ny):
        if height + yp_ic[j] - yp_ic[j-1] <= sl_99_ic:
            npsl += 1 
        height += yp_ic[j] - yp_ic[j-1]
    	
    #!---------------------------------------------!	
    			 
    # Calculation of the number of mesh nodes in the viscous sublayer at cf peak
    npvis = 0     # number of points viscous sublayer
    height = 0.0  # cumulative height in viscous unit (y+)
    
    # Rescaling y-coordinates with viscous unit at peak cf
    yp_peak = yp / delta_nu_peak
      
    for j in range(1, ny):
        if height + yp_peak[j] - yp_peak[j-1] <= 5.0:
            npvis += 1 
        height += yp_peak[j] - yp_peak[j-1]
    
    #!---------------------------------------------!
        
    # Delta y+ at the BL thickness (d: delta) at Re_tau = 500
    c = 0         # integer to store the index (c: counter)
    
    for j in range(1,ny):
        if yp[j] < bl_thickness:
            c = c + 1
    
    # We consider the element just above the estimated BL edge
    delta_yd = yp[c+1] - yp[c]
    
    # Mesh sizes at Re_tau = 500
    delta_x_nd_500  = delta_x  / delta_nu_500
    delta_y1_nd_500 = delta_y1 / delta_nu_500
    delta_z_nd_500  = delta_z  / delta_nu_500
    
    # BL edge
    delta_yd_nd_500 = delta_yd / delta_nu_500
    
    #!---------------------------------------------!  
     
    # Printing useful information to the screen
    print()
    print('!----- Inputs: -----!')
    print()
    print('Number of mesh nodes in streamwise direction,  nx = ', nx)
    print('Number of mesh nodes in wall normal direction, ny = ', ny)
    print('Number of mesh nodes in spanwise direction,    nz = ', nz)
    print()
    print('Domain dimension, Lx = ',xlx)     
    print('Domain dimension, Ly = ',yly)
    print('Domain dimension, Lz = ',zlz)
    print()
    print('Beta parameter = ', beta)
    print('Kinematic viscosity, nu = ', nu)
    print('Wall velocity, Uwall = ', uwall)
    print('Time step, delta_t = ', delta_t)
    print('Tripping wire diameter, twd = ', twd)
    print('Reynolds number (Re = 1/nu) = ', re)
    print()
    print('!--- Reference data according to Cimarelli et al. (2024): ---!')
    print()
    print('Boundary layer thickness at Re_tau = 500, bl_thickness = ', bl_thickness)
    print('Skin friction at peak, cf = ', cf)
    print()
    print('!----- Outputs: -----!')
    print()
    print('!--- Non-dimensional domain dimensions: ---!')
    print()
    print('Length of the domain (Lx+) at IC:', xlx_nd_ic)
    print('Height of the domain (Ly+) at IC:', yly_nd_ic)
    print('Width  of the domain (Lz+) at IC:', zlz_nd_ic)
    print()
    print('Length of the domain (Lx+) at peak cf:', xlx_nd_peak)
    print('Height of the domain (Ly+) at peak cf:', yly_nd_peak)
    print('Width  of the domain (Lz+) at peak cf:', zlz_nd_peak)
    print()
    print('!--- Numerics-related parameters: ---!')
    print()
    print('Estimated CFL at IC:', CFL)
    print('Estimated D   at IC:', D)
    print('Estimated Pé  at IC:', Pe)
    print('Estimated stability parameter S at IC:', S)
    print()
    print()
    print('!--- Mesh sizes at peak cf (cf_max ~ 0.007, Cimarelli et al. (2024)) ---!)
    print()
    print('Mesh size x-direction: delta_x+ =', delta_x_nd_peak)
    print('Mesh size y-direction at the first element near the wall: delta_y1+ =', delta_y1_nd_peak)
    print('Mesh size z-direction: delta_z+ =', delta_z_nd_peak)
    print()
    print('!--- Mesh sizes at Re_tau = 500 (imarelli et al. (2024)) ---!)
    print()
    print('Mesh size x-direction: delta_x+ =', delta_x_nd_500)
    print('Mesh size y-direction at the first element near the wall: delta_y1+ =', delta_y1_nd_500)
    print('Mesh size z-direction: delta_z+ =', delta_z_nd_500)
    print()
    print('Mesh size (y) at the estimated BL edge: delta_yd+ =', delta_yd_nd_500)
    print()
    print('!--- Number of discretization points ---!')
    print()
    
    # to be completed
    print('Number of mesh nodes in the viscous sublayer at cf peak:', npvis)
    print('Number of mesh nodes in the initial shear layer:', npsl)
    print()
    #print('Estimated  initial momentum thickness of the shear layer (approx. 54*nu/U_wall) (dimensional): theta_sl =', theta_sl)
    #print('Calculated initial thickness of the shear layer (y+ where Umean < 0.01 Uwall) (non-dimensional): sl_99^+_IC =', sl_99_ic)
            
    # Create data arrays with inputs
    data = [
            ["nx/ny/nz", "Lx/Ly/Lz" ],
            [nx,          xlx       ],
            [ny,          yly       ],
            [nz,          zlz       ]
           ] 
    
    data2 = [
             ["beta", "cf", "nu",  "Re", "uwall", "delta_t", "twd"],
             [ beta,   cf,   nu,    re,   uwall,   delta_t,   twd ],
            ]

    # Create the tables using tabulate
    table  = tabulate(data,  headers="firstrow", tablefmt="fancy_grid")
    table2 = tabulate(data2, headers="firstrow", tablefmt="fancy_grid")
    
    # Save the tables as a text file 
    with open("sim_settings.txt", "w") as f:
         f.write("!--- Temporal TBL setting parameters ---!\n")
         f.write("\n")
         f.write("!--- Inputs: ---!\n")
         f.write(table)
         f.write("\n")
         f.write(table2)
         f.write("\n")
                 
    # Create data array with outputs
    data = [
            ["Lx+/Ly+/Lz+ at peak cf", "Lx+/Ly+/Lz+ at IC", "CFL/Pé/D", "delta_y1+/delta_yd+/delta_yn+", "delta_x+/delta_z+", "AR_x1/AR_xn", "AR_z1/AR_zn"],
            [ xlx_nd,                   xlx_nd_ic,           CFL,        delta_y1,                        delta_x,             AR_x1,         AR_z1,      ],
            [ yly_nd,                   yly_nd_ic,           Pe,         delta_yd,                        delta_z,             AR_xn,         AR_zn,      ],
            [ zlz_nd,                   zlz_nd_ic,           D,          delta_yn,                        "/",                 "/",           "/"         ],
            ["---",                     "---",               "---",      "---",                           "---",               "---",         "---"       ],
            ["npvis",                   "npsl",              "theta_sl", "sl_99^+_IC",                    "sh_vel_IC",         "sh_vel_peak", "n_tot"     ],
            [ npvis,                     npsl,                theta_sl,   sl_99_ic,                        sh_vel_ic,           sh_vel,        ntot       ],
           ] 

    # Create the table using tabulate
    table = tabulate(data, headers="firstrow", tablefmt="fancy_grid")

    # Save the table as a text file 
    with open("sim_settings.txt", "a") as f:
         f.write("\n")
         f.write("!--- Outputs: ---!\n")
         f.write(table)
            
    # Opening to add final informations
    with open('sim_settings.txt', 'a') as f: 
         f.write("\n")
         f.write("\n")
         f.write("!--- INFO: ---!\n")
         f.write("We are employing 2 different cf values:\n")
         f.write("\n")
         f.write("1) Value obtained from the derivative at the wall due to the IC,\n")
         f.write("   and it is used in order to estimate the number of points in the initial shear layer and the domain dimensions at IC.\n")
         f.write("\n")     
         f.write("2) Peak value found in literature (e.g. 0.007, see Cimarelli et al. (2024)),\n")
         f.write("   and it is used to evaluate the spacings of the mesh (delta_x+,delta_y+,delta_z+).\n")
         f.write("\n")
         f.write("Stability parameter, S < 1: "+ str(S) +" (Thompson et al. (1985))\n")
         f.write("\n")
         f.write("!--- List of acronyms & variables: ---!\n")
         f.write("AR:            Aspect Ratio of mesh elements (width/height).\n")
         f.write("npvis:         Number of points viscous sublayer at cf peak (y+ < 5).\n")
         f.write("npsl:          Number of points initial shear layer (y+ < theta_sl_true+).\n")
         f.write("theta_sl:      Estimated  initial momentum thickness of the shear layer (approx. 54*nu/U_wall) (dimensional).\n")
         f.write("sl_99^+_IC:    Calculated initial thickness of the shear layer (y+ where Umean < 0.01 Uwall) (non-dimensional).\n")
         f.write("sh_vel_IC:     Shear velocity of the initial condition.\n")
         f.write("sh_vel_peak:   Shear velocity at peak cf, according to Cimarelli et al. (2024).\n")
         f.write("!-------------------------------------!\n")
         
         
    


    
    



