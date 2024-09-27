
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: This file stores Python subroutines for the mesh of 
!              Incompact3d, used in the setup of simulations.
!              Subroutines list:
!               - stretching_mesh_y;
!               - calculate_geometric_quantities. 
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
config_path = os.path.abspath(os.path.join(current_dir, '../', '3-post_processing/python_common'))
sys.path.append(config_path)

# Import the plotting_params module
import plot_params as pp

#!-----------------------------------------------------------------------------!

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Adaptation of original 'stretching' subroutine of Incompact3d. 
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def stretching_mesh_y(ny, yly, beta, istret):

    # Declare local constant Pi
    pi = np.pi

    # Work variables
    yeta = np.zeros(ny)
    yp   = np.zeros(ny)
    
    # If periodic BC is imposed along y, nym = ny, otherwise nym = ny - 1
    nym = ny - 1 
    
    #!-----------------------------------------------------------------------------!
    
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
    
    return yp

    #!-----------------------------------------------------------------------------!

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculate and write to .txt file geometric quantities: 
!               - y-coordinates of elements' faces (yp);
!               - height of elements in y-direction (delta_y);
!               - Growth-Rate in y-direction (GR_y);
!               - Aspect-Ratio in xy plane (AR_xy).
!              Used to check yp values with default ones of Incompact3d and
!              to further analyse the geometric properties of the elements.
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def calculate_geometric_quantities(ny, yp, delta_x):

    # Delta y
    delta_y = np.zeros(ny)

    # Growth-Rate (GR) of grid elements in y-direction
    gr_y = np.zeros(ny)

    # AR xy-plane
    AR_xy = np.zeros(ny)
    
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
    
    #!--- Writing the results to .txt files ---!

    # Write yp, delta_y and GR_y in a .txt file
    # (checked, yp values are equal to the ones printed by Incompact3d)

    with open('mesh_y.txt', 'w') as f:
        f.write(f"{'yp':>{pp.c_w}}, "      +
                f"{'delta_y':>{pp.c_w}}, " +
                f"{'gr_y':>{pp.c_w}}, "    +
                f"{'AR_xy':>{pp.c_w}}\n"   )
            
        for j in range(0,ny):
            f.write(f"{yp[j]:{pp.fs6}}, "      +
                    f"{delta_y[j]:{pp.fs6}}, " + 
                    f"{gr_y[j]:{pp.fs6}}, "    +
                    f"{AR_xy[j]:{pp.fs6}}\n"   )

    #!-----------------------------------------------------------------------------!


