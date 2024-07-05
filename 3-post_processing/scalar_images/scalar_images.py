#!---------------------------------------------------------!
#! With this script, we perform plotting of 2D z-normal    !
#! planes of instantaneous scalar field, used for visual   !
#! check of the evolution of a TTBL.                       !
#!---------------------------------------------------------!

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '..', 'python_common'))
sys.path.append(config_path)

# Import the plotting_params module
import plot_params as pp

# Import function to set plots
from plot_settings import set_plot_settings

# Import function to read 'input.i3d' and 'post.prm' files
from read_incompact3d_files import read_input_files

#!--------------------------------------------------------------------------------------!

# Create folders to store images
os.makedirs('images', mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
itype, nx, ny, nz, Lx, Ly, Lz, re, iswitch_wo, file1, filen, icrfile, nr, add_string = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!

# Create x-coordinates vector
x = np.linspace(0.0, Lx, nx)

# Read y-coordinates vector
y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)

# Define the binary file path
file_path = 'data/planes/phiplanez-0059.bin'


# Read the binary file into a numpy array
with open(file_path, 'rb') as file:
    data = np.fromfile(file, dtype=np.float64)

# Reshape the fields to 2D arrays using Fortran order
data = data.reshape((nx, ny), order='F')


# Create meshgrid from x and y
X, Y = np.meshgrid(x, y)

# Mean velocity profile
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

# Limits for axes
xliminf = 0.0
xlimsup = Lx
yliminf = 0.0
ylimsup = Ly

# Plotting with contourf
ax.contourf(X, Y, data.T, cmap='viridis')  # Plotting directly with X, Y, and data
ax.colorbar(label='Value')  # Add a colorbar with label

# Axes labels
ax.set_xlabel(r'$x/D$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$y/D$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)

# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Saving the figure
plt.savefig(f'images/phiplanez-{snap_numb}_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)
    
# Show the figure
plt.show()





