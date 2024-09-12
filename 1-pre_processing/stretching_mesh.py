"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Adaptation of original 'stretching' subroutine of Incompact3d. 
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

import numpy as np

def stretching_mesh_y(ny, yly, beta, istret):

    # Declare local constant Pi
    pi = np.pi

    # Work variables
    yeta = np.zeros(ny)
    yp   = np.zeros(ny)

    # If periodic BC is imposed along y, nym = ny, otherwise nym = ny - 1
    nym = ny - 1 

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

