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

# Print the reshaped data (optional)
print(data)


# Plotting with contourf
plt.figure(figsize=(8, 6))
plt.contourf(x, y, data, cmap='viridis')  # Plotting directly with X, Y, and data
plt.colorbar(label='Value')  # Add a colorbar with label
plt.title('Contour Plot of Data')
plt.xlabel('X Axis')
plt.ylabel('Y Axis')
plt.grid(True)
plt.show()





