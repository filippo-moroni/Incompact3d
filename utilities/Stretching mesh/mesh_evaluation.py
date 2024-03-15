#!---------------------------------------------------------!
#! This is the Python version of 'mesh_evaluation.f90'     !
#! Same values are obtained from both codes (tested)       !
#!                                                         !
#! Y-coordinates at the center of mesh nodes are not       !
#! evaluated here.                                         !
#!---------------------------------------------------------!

import numpy as np

# Inputs
xlx = 32.0   # domain dimension in x direction
nx = 64      # number of points in x direction
yly = 8.0    # domain dimension in y direction
ny = 129     # number of points in y direction
zlz = 16.0   # domain dimension in z direction
nz = 64      # number of points in z direction

nym = ny - 1 # if periodic BC is imposed, nym = ny, otherwise nym = ny - 1

istret = 3   # y mesh refinement (0:no, 1:center, 2:both sides, 3:bottom)
beta = 0.4   # beta parameter for mesh stretching
cf = 0.007   # maximum cf estimated
nu = 0.002   # kinematic viscosity (if D = 1 and U_wall = 1, Re_D = 500)
uwall = 1.0  # velocity of the wall

# Declare local constant Pi
pi = np.pi

# Work variables
yeta = np.zeros(ny)
yp = np.zeros(ny)

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
    
    # Rescaling also domain dimensions
    xlx_nd = xlx / delta_nu
    zlz_nd = zlz / delta_nu

    # First and last elements' dimension
    delta_y1 = yp[2] - yp[1]
    delta_yn = yp[ny-1] - yp[ny-2]
    
    # Calculation of the number of mesh nodes in the viscous sublayer
    npvis = 0     # number of points viscous sublayer
    height = 0.0  # cumulative height in viscous unit (y+)
  
    for j in range(1, ny):
        if height + yp[j] - yp[j-1] <= 5:
            npvis += 1 
        height += yp[j] - yp[j-1]

    # Printing useful information to the screen
    print()
    print('!--- Inputs: ---!')
    print()
    print('Domain dimension, Lx = ',xlx)
    print('Number of mesh nodes in streamwise direction, nx = ', nx)  
    print('Domain dimension, Ly = ',yly)
    print('Number of mesh nodes in wall normal direction, ny = ', ny)
    print('Domain dimension, Lz = ',zlz)
    print('Number of mesh nodes in spanwise direction, nz = ', nz)
    print()
    print('Beta parameter = ', beta)
    print('Skin friction coefficient employed, cf = ', cf)
    print('Kinematic viscosity, nu = ', nu)
    print('Wall velocity, Uwall = ', uwall)
    print()
    print('!--- Outputs: ---!')
    print('Mesh size at the first element near the wall: delta_y1+ =', delta_y1)
    print('Mesh size at the last element away from the wall: delta_yn+ =', delta_yn)
    print('Number of mesh nodes in the viscous sublayer:', npvis)
    print()
    print('Length of the domain (Lx+):', xlx_nd)
    print('Height of the domain (Ly+):', yp[ny-1])
    print('Width of the domain (Lz+):', zlz_nd)
    print()
    print('Mesh size x-direction: delta_x+ =', delta_x)
    print('Mesh size z-direction: delta_z+ =', delta_z)
    
    
    



