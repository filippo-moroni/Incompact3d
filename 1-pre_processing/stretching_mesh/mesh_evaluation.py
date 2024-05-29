#!---------------------------------------------------------!
#! mesh_evaluation.py:                                     !
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
# flow parameters for temporal TBL simulations and for
# channel flow simulations
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
istret = 2               # y mesh refinement (0:no, 1:center, 2:both sides, 3:bottom)
beta = 0.25              # beta parameter for mesh stretching
delta_t = 0.002          # time-step

# Reynolds number
re = 4764.0              # TTBL: if D = 1 and U_wall = 1, re = Re_D = 500; Channel: re = Re_0 of a laminar Poiseuille flow
nu = 1.0/re              # kinematic viscosity as defined in Incompact3d

# Friction Reynolds number for a channel, considering the centerline Reynolds number of a laminar Poiseuille flow
re_tau = 0.116*re**0.88

# TTBL inputs
twd = 1.0                # trip wire diameter D
bl_thickness = 21.7*twd  # displacement thickness of the temporal TBL at Re_tau = 500 (Cimarelli et al. (2024))

# Reference velocity
uref = 0.667             # reference velocity (TTBL: wall velocity; Channel: bulk velocity)
uwall = uref

# Selection of friction coefficient at peak or at steady state conditions
# for grid spacings evaluation in wall units 

cf = 0.007               # maximum cf estimated at peak (Cimarelli et al. (2024))
cf = 0.00793             # steady state cf of a channel flow at Re_tau = 200 (Quadrio & Ricco (2004)) 

nx = 36                  # number of points in x direction
ny = 649                 # number of points in y direction
nz = 54                  # number of points in z direction

xlx = 10.0               # domain dimension in x direction
#yly = 3.0*bl_thickness   # domain dimension in y direction
yly = 2.0                # domain dimension in y direction
zlz = 8.0                # domain dimension in z direction

nym = ny - 1             # if periodic BC is imposed along y, nym = ny, otherwise nym = ny - 1

# Declare local constant Pi
pi = np.pi

# Work variables
yeta = np.zeros(ny)
yp   = np.zeros(ny)
Uo   = np.zeros(ny)

# Delta y
delta_y = np.zeros(ny)

# Growth-Rate (GR) of grid elements in y-direction
gr_y = np.zeros(ny)

# AR xy-plane
AR_xy = np.zeros(ny)

# Total number of points
n_tot = nx*ny*nz

# Start of calculations
yinf = -yly/2.0
den = 2.0 * beta * yinf
xnum = -yinf - np.sqrt(pi*pi*beta*beta + yinf*yinf)
alpha = np.abs(xnum/den)
xcx = 1.0/beta/alpha

# First point initialization
if istret == 1:
    yp[0] = 0.0
if istret == 2:
    yp[0] = 0.0
if istret == 1:
    yeta[0] = 0.0
if istret == 2:
    yeta[0] = -0.5
if istret == 3:
    yp[0] = 0.0
if istret == 3:
    yeta[0] = -0.5

# Stretched mesh
if alpha != 0.0:
    # Loop from j = 1 to ny - 1
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

# Uniform mesh
if alpha == 0.0:
    yp[0] = -1.e10
    for j in range(1, ny):
        yeta[j] = j*(1.0/ny)
        yp[j] = -beta * np.cos(pi*yeta[j]) / np.sin(yeta[j]*pi)


#!--- Calculations valid for both TTBL and Channel ---!

# Calculate the spacings along x and z (uniform)
delta_x = xlx / nx
delta_z = zlz / nz

# Shear velocity at peak cf or at steady state 
sh_vel_peak = np.sqrt((cf/2.0)) * uref

# Viscous length at peak cf or at steady state
delta_nu_peak = nu / sh_vel_peak

# First element y-dimension
delta_y1 = yp[1] - yp[0]

# Rescaling the mesh spacings with viscous unit at peak cf or at steady state
delta_x_nd_peak  = delta_x  / delta_nu_peak
delta_y1_nd_peak = delta_y1 / delta_nu_peak
delta_z_nd_peak  = delta_z  / delta_nu_peak

# Non-dimensional domain dimensions at peak cf or at steady state (nd: non-dimensional)
xlx_nd_peak = xlx / delta_nu_peak
yly_nd_peak = yly / delta_nu_peak
zlz_nd_peak = zlz / delta_nu_peak

#!--- Estimation of numerics-related parameters (CFL, D, Pé, S) at IC or at steady state ---!
   
CFL = uref * delta_t / delta_x      
D =   nu * delta_t / (delta_y1**2)    
Pe =  uref * delta_x / nu
        
# Stability parameter (S < 1) (see Thompson et al. (1985)) 
S = ((uref**2)*delta_t)/(2.0*nu)

#!--- Calculation of y-dir. geometric quantities ---!

# Delta y of grid elements
for j in range(1,ny):
    delta_y[j] = yp[j] - yp[j-1]
            
# Calculate the GR (Growth-Rate) of the elements in y-direction
for j in range(2,ny):
    gr_y[j] = delta_y[j]/delta_y[j-1]
        
# AR of the elements in xy plane (width / heigth)
for j in range(1,ny):
    AR_xy[j] = delta_x / delta_y[j]

#!--------------------------------------------------!

# This part is valid for meshes with refinement at the bottom boundary and for temporal TBL cases 
if istret == 3:
   
    # Calculating the initial thickness of the shear layer (see Kozul et al. (2016))
    theta_sl = 54.0 * nu / uwall
    
    # Mean gradient due to initial condition (analytical derivative)  
    mg = - uwall / (4.0 * theta_sl) * (1.0 / np.cosh(twd / 2.0 / theta_sl))**2
    
    #!--- Extra shear velocities for TTBL (IC and Re_tau = 500): ---!
    
    # Shear velocity due to initial condition
    sh_vel_ic = np.sqrt(nu * np.abs(mg))
       
    # Shear velocity at Re_tau = 500:
    sh_vel_500 = 500.0 * nu / bl_thickness
    
    #!--------------------------------------------------------!

    # Calculate viscous length at IC and at Re_tau = 500
    delta_nu_ic   = nu / sh_vel_ic
    delta_nu_500  = nu / sh_vel_500

    # Rescaling the mesh spacings with viscous unit at IC
    delta_x_nd_ic  = delta_x  / delta_nu_ic
    delta_y1_nd_ic = delta_y1 / delta_nu_ic
    delta_z_nd_ic  = delta_z  / delta_nu_ic
            
    # Non-dimensional domain dimensions at IC (nd: non-dimensional)  
    xlx_nd_ic = xlx / delta_nu_ic
    yly_nd_ic = yly / delta_nu_ic
    zlz_nd_ic = zlz / delta_nu_ic
       
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
    
    # We consider the element just below the estimated BL edge
    delta_yd = yp[c] - yp[c-1]
    
    # Mesh sizes at Re_tau = 500
    delta_x_nd_500  = delta_x  / delta_nu_500
    delta_y1_nd_500 = delta_y1 / delta_nu_500
    delta_z_nd_500  = delta_z  / delta_nu_500
    
    # BL edge
    delta_yd_nd_500 = delta_yd / delta_nu_500
    
    #!---------------------------------------------!
             
#!--- Printing useful information to the screen ---!
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
print('Stretching index =', istret)
print('Beta parameter = ', beta)
print('Kinematic viscosity, nu = ', nu)
print('Time step, delta_t = ', delta_t)
print('Reynolds number (Re = 1/nu) = ', re)

if istret != 3:
    print('Reference velocity U_ref = ', uref)
    print()
    print('Skin friction coefficient at steady state, cf = ', cf)
    print()
    print('Estimated friction Reynolds number, Re_tau ~ ', re_tau)
    print()
    
elif istret == 3:
    print('Wall velocity, Uwall = ', uwall)
    print('Tripping wire diameter, twd = ', twd)
    print('!--- Reference data according to Cimarelli et al. (2024): ---!')
    print()
    print('Boundary layer thickness at Re_tau = 500, bl_thickness = ', bl_thickness)
    print('Skin friction coefficient at peak, cf = ', cf)
    
print()
print()
print('!----- Outputs: -----!')
print()
print('!--- Non-dimensional domain dimensions: ---!')
print()

if istret == 3:
    print('Length of the domain (Lx+) at IC:', xlx_nd_ic)
    print('Height of the domain (Ly+) at IC:', yly_nd_ic)
    print('Width  of the domain (Lz+) at IC:', zlz_nd_ic)
    print()

print('Length of the domain (Lx+) at peak cf or at steady state:', xlx_nd_peak)
print('Height of the domain (Ly+) at peak cf or at steady state:', yly_nd_peak)
print('Width  of the domain (Lz+) at peak cf or at steady state:', zlz_nd_peak)
print()
print('!--- Numerics-related parameters based on reference velocity U_ref: ---!')
print()
print('Estimated CFL,x:', CFL)
print('Estimated D,y  :', D)
print('Estimated Pé,x :', Pe)
print('Estimated stability parameter S,x:', S)
print()
print('!--- Mesh sizes at peak cf or at steady state ---!')
print()
print('Mesh size x-direction: delta_x+ =', delta_x_nd_peak)
print('Mesh size y-direction at the first element near the wall: delta_y1+ =', delta_y1_nd_peak)
print('Mesh size z-direction: delta_z+ =', delta_z_nd_peak)
print()

if istret == 3:
    print('!--- Mesh sizes at Re_tau = 500 (Cimarelli et al. (2024)) ---!')
    print()
    print('Mesh size x-direction: delta_x+ =', delta_x_nd_500)
    print('Mesh size y-direction at the first element near the wall: delta_y1+ =', delta_y1_nd_500)
    print('Mesh size z-direction: delta_z+ =', delta_z_nd_500)
    print()
    print('Mesh size (y) at the estimated BL edge: delta_yd+ =', delta_yd_nd_500)
    print()
    print('!--- Number of discretization points ---!')
    print()
    print('Number of mesh nodes in the viscous sublayer at cf peak:', npvis)
    print('Number of mesh nodes in the initial shear layer:', npsl)
    print()
    
    #print('Estimated  initial momentum thickness of the shear layer (approx. 54*nu/U_wall) (dimensional): theta_sl =', theta_sl)
    #print('Calculated initial thickness of the shear layer (y+ where Umean < 0.01 Uwall) (non-dimensional): sl_99^+_IC =', sl_99_ic)
    
#!-------------------------------------------------!

# Writing the results also into files
  
# Column width for writing to .txt file
c_w = 24  

# Write yp, delta_y and GR_y in a .txt file
with open('mesh_y.txt', 'w') as f:
    f.write(f"{'yp':<{c_w}}, {'delta_y':<{c_w}}, {'gr_y':<{c_w}}, {'AR_xy':<{c_w}}\n")
    for j in range(0,ny):
        f.write(f"{yp[j]:<{c_w}}, {delta_y[j]:<{c_w}}, {gr_y[j]:<{c_w}}, {AR_xy[j]:<{c_w}}\n")
     
# Create data arrays with inputs
data = [
        ["nx/ny/nz", "Lx/Ly/Lz" ],
        [nx,          xlx       ],
        [ny,          yly       ],
        [nz,          zlz       ]
       ] 

data2 = [
         ["beta", "nu", "U_ref", "delta_t", "Re", "cf"],
         [ beta,   nu,   uref,    delta_t,   re,   cf ],
        ]

# TTBL file creation and saving
if istret == 3:    
    data2 = [
             ["beta", "nu", "uwall", "delta_t", "twd", "Re"],
             [ beta,   nu,   uwall,   delta_t,   twd,   re ],
            ]
    data3 = [
             ["bl_thickness", "cf_max"],
             [ bl_thickness,   cf     ],
            ]
                   
    # Create the tables using tabulate
    table  = tabulate(data,  headers="firstrow", tablefmt="fancy_grid")
    table2 = tabulate(data2, headers="firstrow", tablefmt="fancy_grid")
    table3 = tabulate(data3, headers="firstrow", tablefmt="fancy_grid")
    
    # Save the tables as a text file 
    with open("sim_settings.txt", "w") as f:
         f.write("!----- Temporal TBL setting parameters -----!\n")
         f.write("\n")
         f.write("!----- Inputs: -----!\n")
         f.write(table)
         f.write("\n")
         f.write(table2)
         f.write("\n")
         f.write("!--- BL thickness @ Re_tau = 500 and cf peak, according to Cimarelli et al. (2024) ---!\n")
         f.write(table3)
                 
    # Create data arrays with outputs
    data = [
            ["Lx+/Ly+/Lz+ at IC", "Lx+/Ly+/Lz+ at peak cf" ],
            [ xlx_nd_ic,           xlx_nd_peak             ],
            [ yly_nd_ic,           yly_nd_peak             ],
            [ zlz_nd_ic,           zlz_nd_peak             ],
           ]
           
    data2 = [
             ["CFL,x", "D,y", "Pé,x", "S,x"],
             [ CFL,     D,     Pe,     S   ],
            ]
            
    data3 = [
             ["/",           "IC",           "peak cf",        "Re_tau = 500"   ],
             ["delta_x+",     delta_x_nd_ic,  delta_x_nd_peak,  delta_x_nd_500  ],
             ["delta_y1+",    delta_y1_nd_ic, delta_y1_nd_peak, delta_y1_nd_500 ],
             ["delta_z+",     delta_z_nd_ic,  delta_z_nd_peak,  delta_z_nd_500  ],
             ["delta_yd+",   "/",            "/",               delta_yd_nd_500 ],       
            ]
           
    data4 = [
             ["npvis", "npsl", "theta_sl", "sl_99^+_IC", "sh_vel_IC", "sh_vel_peak", "sh_vel Re_tau = 500", "n_tot" ],
             [ npvis,   npsl,   theta_sl,   sl_99_ic,     sh_vel_ic,   sh_vel_peak,   sh_vel_500,            n_tot  ],
            ] 

    # Create the tables using tabulate
    table  = tabulate(data,  headers="firstrow", tablefmt="fancy_grid")
    table2 = tabulate(data2, headers="firstrow", tablefmt="fancy_grid")
    table3 = tabulate(data3, headers="firstrow", tablefmt="fancy_grid")
    table4 = tabulate(data4, headers="firstrow", tablefmt="fancy_grid")

    # Save the table as a text file and final informations
    with open("sim_settings.txt", "a") as f:
         f.write("\n")
         f.write("!----- Outputs: -----!\n")
         f.write("\n")
         f.write("!--- Non-dimensional domain dimensions: ---!\n")
         f.write(table) 
         f.write("\n")
         f.write("!--- Numerics-related parameters: ---!\n")
         f.write(table2) 
         f.write("\n")
         f.write("!--- Mesh sizes at IC, at peak cf and at Re_tau = 500: ---!\n")
         f.write(table3) 
         f.write("\n")
         f.write("!--- Miscellaneous ---!\n")
         f.write(table4) 
         f.write("\n")
         f.write("\n")
         f.write("!--- INFO: ---!\n")
         f.write("We are employing 3 different cf values:\n")
         f.write("\n")
         f.write("1) Value obtained from the derivative at the wall due to the IC.\n")
         f.write("\n")     
         f.write("2) Peak value found in literature (e.g. 0.007, see Cimarelli et al. (2024)).\n")
         f.write("\n")
         f.write("3) Value at Re_tau = 500, again according to Cimarelli et al. (2024).\n")
         f.write("\n")
         f.write("!--- List of acronyms & variables: ---!\n")
         f.write("\n")
         f.write("S:             Stability parameter, S < 1 (Thompson et al. (1985)).\n")
         f.write("npvis:         Number of points viscous sublayer at cf peak (y+ < 5).\n")
         f.write("npsl:          Number of points initial shear layer (y+ < sl_99^+_IC).\n")
         f.write("theta_sl:      Estimated  initial momentum thickness of the shear layer (approx. 54*nu/U_wall) (dimensional).\n")
         f.write("sl_99^+_IC:    Calculated initial thickness of the shear layer (y+ where Umean < 0.01 Uwall) (non-dimensional).\n")
         f.write("sh_vel_IC:     Shear velocity of the initial condition.\n")
         f.write("sh_vel_peak:   Shear velocity at peak cf, according to Cimarelli et al. (2024).\n")
         f.write("sh_vel_500:    Shear velocity at Re_tau = 500, according to Cimarelli et al. (2024).\n")
         f.write("\n")
         f.write("!-------------------------------------!\n")
         
# Channel file creation and saving         
else:

    # Create the tables using tabulate
    table  = tabulate(data,  headers="firstrow", tablefmt="fancy_grid")
    table2 = tabulate(data2, headers="firstrow", tablefmt="fancy_grid")
    
    # Save the tables as a text file 
    with open("sim_settings.txt", "w") as f:
         f.write("!----- Channel setting parameters -----!\n")
         f.write("\n")
         f.write("!----- Inputs: -----!\n")
         f.write(table)
         f.write("\n")
         f.write(table2)
    
    # Create data arrays with outputs
    data = [
            ["Lx+/Ly+/Lz+", "dx+/dy+/dz+"     ],
            [ xlx_nd_peak,   delta_x_nd_peak  ],
            [ yly_nd_peak,   delta_y1_nd_peak ],
            [ zlz_nd_peak,   delta_z_nd_peak  ],
           ]
           
    data2 = [
             ["CFL,x", "D,y", "Pé,x", "S,x"],
             [ CFL,     D,     Pe,     S   ],
            ]
                       
    data3 = [
             ["sh_vel",    "n_tot" ],
             [ sh_vel_peak, n_tot  ],
            ]
            
    # Create the tables using tabulate
    table  = tabulate(data,  headers="firstrow", tablefmt="fancy_grid")
    table2 = tabulate(data2, headers="firstrow", tablefmt="fancy_grid")
    table3 = tabulate(data3, headers="firstrow", tablefmt="fancy_grid")
    
    # Save the table as a text file and final informations
    with open("sim_settings.txt", "a") as f:
         f.write("\n")
         f.write("!----- Outputs: -----!\n")
         f.write("\n")
         f.write("!--- Non-dimensional domain dimensions and grid spacings: ---!\n")
         f.write(table) 
         f.write("\n")
         f.write("!--- Numerics-related parameters: ---!\n")
         f.write(table2) 
         f.write("\n")
         f.write("!--- Shear velocity and total number of points: ---!\n")
         f.write(table3) 
         f.write("\n")
         f.write("!-------------------------------------!\n") 

         
         
    


    
    



