#!---------------------------------------------------------!
#! With this script, we perform plotting of 2D z-normal    !
#! planes of instantaneous scalar field, used for visual   !
#! check of the evolution of a TTBL.                       !
#!---------------------------------------------------------!

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

#!--- Mesh section ---!

# Create x-coordinates vector
x = np.linspace(0.0, Lx, nx)

# Read y-coordinates vector
y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)

# Create meshgrid from x and y
X, Y = np.meshgrid(x, y)

#!--------------------------------------------------------------------------------------!

# Define the binary file path
file_path = 'data/planes/phiplanez-0059.bin'


# Read the instantaneous scalar field binary file into a numpy array
with open(file_path, 'rb') as file:
    data = np.fromfile(file, dtype=np.float64)

# Reshape scalar field to 2D array using Fortran order
data = data.reshape((nx, ny), order='F')

# Transpose and flip upside-down
data = np.flipud(data.T)

# Limits for axes (not used in reality, but needed to use 'set_plot_settings'
xliminf = 0.0
xlimsup = 1.0
yliminf = 0.0
ylimsup = 1.0

# Subplots environment
fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)


# Extent of the image (dimensions of the domain)
extent = [0.0, Lx, 0.0, Ly]

# Plotting 
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='10%', pad=0.05)

im = ax.imshow(data, cmap='Blues', extent=extent, vmin=0.0, vmax=1.0)

fig.colorbar(im, cax=cax, orientation='vertical')

# Axes labels
ax.set_xlabel(r'$x/D$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
ax.set_ylabel(r'$y/D$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)


# Set the plot parameters using the function 'set_plot_settings'
# Last argument is the switcher for semilog plot (1: yes, 0: no)
set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

# Saving the figure
plt.savefig(f'images/phiplanez_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)
    
# Show the figure
plt.show()





