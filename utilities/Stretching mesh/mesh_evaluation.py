#!---------------------------------------------------------!
#! This an improved Python version of 'mesh_evaluation.f90'!
#! Same values are obtained from both codes (tested)       !
#!                                                         !
#! Y-coordinates at the center of mesh nodes are not       !
#! evaluated here.                                         !
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
nx = 64           # number of points in x direction
ny = 129          # number of points in y direction
nz = 64           # number of points in z direction

xlx = 8.0         # domain dimension in x direction
yly = 4.0         # domain dimension in y direction
zlz = 8.0         # domain dimension in z direction

nym = ny - 1      # if periodic BC is imposed, nym = ny, otherwise nym = ny - 1

istret = 3        # y mesh refinement (0:no, 1:center, 2:both sides, 3:bottom)
beta = 0.3        # beta parameter for mesh stretching
cf = 0.007        # maximum cf estimated
nu = 0.002        # kinematic viscosity (if D = 1 and U_wall = 1, Re_D = 500)
uwall = 1.0       # velocity of the wall
delta_t = 0.05    # time-step
twd = 1.0         # trip wire diameter D

re = 1//nu        # Reynolds number as defined in Incompact3d

# Declare local constant Pi
pi = np.pi

# Work variables
yeta = np.zeros(ny)
yp = np.zeros(ny)
Uo = np.zeros(ny)

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

# This part is valid for meshes with refinement at the bottom boundary only
if istret == 3:
    
    # Calculate the spacings along x and z (uniform)
    delta_x = xlx / nx
    delta_z = zlz / nz
    
    # Calculate shear velocity (temporal BL)
    sh_vel = np.sqrt((cf/2.0)) * uwall

    # Calculate viscous unit
    delta_nu = nu / sh_vel

    # Rescaling the y coordinates and the spacings along x and z
    yp = yp / delta_nu
    delta_x = delta_x / delta_nu
    delta_z = delta_z / delta_nu
    
    # Rescaling also domain dimensions (nd: non-dimensional)
    xlx_nd = xlx / delta_nu
    yly_nd = yly / delta_nu
    zlz_nd = zlz / delta_nu  

    # First and last elements' dimension
    delta_y1 = yp[2] - yp[1]
    delta_yn = yp[ny-1] - yp[ny-2]
    
    # Aspect ratios (ARs) of 1st and last elements
    AR_x1 = delta_x / delta_y1	# AR x direction, 1st element
    AR_xn = delta_x / delta_yn	# AR x direction, nth element
    
    AR_z1 = delta_z / delta_y1	# AR z direction, 1st element
    AR_zn = delta_z / delta_yn	# AR z direction, nth element     
    
    # Estimation of parameters of numerics (CFL, Pe, D) at the beginning of simulation
    CFL = uwall * delta_t / delta_x
    Pe =  uwall * delta_x / nu
    D =   nu * delta_t / (delta_x**2)
    
    # Calculating the initial thickness of the shear layer (see Kozul et al. (2016))
    theta_sl = 54.0 * nu / uwall
    
    # Calculating the initial velocity profile
    for j in range(0, ny):
    	Uo[j] = uwall * (0.5 + 0.5 * (math.tanh(twd/2.0/theta_sl*(1.0 - yp[j]/twd))))
    
    # Rescaling the velocity profile and the shear layer thickness
    Uo = Uo / sh_vel
    theta_sl = theta_sl / delta_nu
    
    # Plotting of the initial velocity profile in wall units
    plt.scatter(yp[0:15], Uo[0:15])
    plt.title("Initial velocity profile near the wall", fontsize=30)
    plt.xlabel("$y^+$", fontsize=30)
    plt.ylabel("$U_o^+$", fontsize=30)
    plt.show()
    
    # Calculate the real value of the initial shear layer
    j = 0
    while Uo[j] > Uo[0]*0.01:
    	theta_sl_true = yp[j]
    	j = j + 1
    		 
    # Calculation of the number of mesh nodes in the viscous sublayer
    npvis = 0     # number of points viscous sublayer
    height = 0.0  # cumulative height in viscous unit (y+)
  
    for j in range(1, ny):
        if height + yp[j] - yp[j-1] <= 5:
            npvis += 1 
        height += yp[j] - yp[j-1]
    
    # Calculation of the number of mesh nodes in the initial shear layer
    npsl = 0      # number of points shear layer
    height = 0.0  # cumulative height in viscous unit (y+)
  
    for j in range(1, ny):
        if height + yp[j] - yp[j-1] <= theta_sl_true:
            npsl += 1 
        height += yp[j] - yp[j-1]
    
     
    # Printing useful information to the screen
    print()
    print('!--- Inputs: ---!')
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
    print('Skin friction coefficient employed, cf = ', cf)
    print('Kinematic viscosity, nu = ', nu)
    print('Reynolds number (Re = 1/nu) = ', re) 
    print('Wall velocity, Uwall = ', uwall)
    print('Time step, delta_t = ', delta_t)
    print()
    print('!--- Outputs: ---!')
    print()
    print('Mesh size at the first element near the wall: delta_y1+ =', delta_y1)
    print('Mesh size at the last element away from the wall: delta_yn+ =', delta_yn)
    print('Number of mesh nodes in the viscous sublayer:', npvis)
    print('Number of mesh nodes in the initial shear layer:', npsl)
    print()
    print('Length of the domain (Lx+):', xlx_nd)
    print('Height of the domain (Ly+):', yly_nd)
    print('Width  of the domain (Lz+):', zlz_nd)
    print()
    print('Mesh size x-direction: delta_x+ =', delta_x)
    print('Mesh size z-direction: delta_z+ =', delta_z)
    print()
    print('Aspect ratio x-direction 1st element: AR_x1 =', AR_x1)
    print('Aspect ratio x-direction nth element: AR_xn =', AR_xn)
    print('Aspect ratio z-direction 1st element: AR_z1 =', AR_z1)
    print('Aspect ratio z-direction nth element: AR_zn =', AR_zn)
    print()
    print('Estimated CFL at t = 0:', CFL)
    print('Estimated Pé  at t = 0:', Pe)
    print('Estimated D   at t = 0:', D)
    print()
    print('Estimated  initial thickness of the shear layer: theta_sl+ =', theta_sl)
    print('Calculated initial thickness of the shear layer: theta_sl_true+ =', theta_sl_true)
    
    # Export inputs and outputs in a .txt file
    data = [
            ["nx/ny/nz", "Lx/Ly/Lz", "Lx+/Ly+/Lz+", "CFL/Pé/D", "AR_x1/AR_xn", "AR_z1/AR_zn", "delta_y1+/delta_yn+", "delta_x+/delta_z+"],
            [nx,          xlx,        xlx_nd,        CFL,        AR_x1,         AR_z1,         delta_y1,              delta_x           ],
            [ny,          yly,        yly_nd,        Pe,         AR_xn,         AR_zn,         delta_yn,              delta_z           ],
            [nz,          zlz,        zlz_nd,        D,          "/",           "/",           "/",                   "/"               ],
            ["---",       "---",      "---",         "---",      "---",         "---",         "---",                 "---"             ],
            ["beta",      "cf",       "nu",          "Re",       "uwall",       "delta_t",     "npvis",               "npsl",           ],
            [ beta,        cf,         nu,            re,         uwall,         delta_t,       npvis,                 npsl,            ],
           ] 

    # Create the table using tabulate
    table = tabulate(data, headers="firstrow", tablefmt="fancy_grid")

    # Save the table as a text file 
    with open("sim_settings.txt", "w") as f:
         f.write(table)
    


    
    



