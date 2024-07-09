#!---------------------------------------------------------!
#! With this script, we perform plotting and saving of     !
#! instantaneous planes of scalar field or of streamwise   !
#! vorticity.                                              !
#!                                                         !
#! Inspired by 'snap2png.py' by R. Corsini                 !
#!---------------------------------------------------------!

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
import matplotlib

# External component to avoid memory leaks in recursive plotting
matplotlib.use('Agg')

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

# Import function to read friction Reynolds number Re_tau from .xdmf files
from read_retau import extract_re_tau_value 

#!--------------------------------------------------------------------------------------!

# Create folders to store images
os.makedirs('images', mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
itype, nx, ny, nz, Lx, Ly, Lz, re, iswitch_wo, file1, filen, icrfile, nr, add_string = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!

# Asking the user what he wants to plot (scalar field, streamwise vorticity)
print()
switcher = int(input("Specify the selector for plotting (0: scalar field, 1: streamwise vorticity): "))
print()

# Scalar field
if switcher == 0:

    Lxi         = Lx
    nxi         = nx
    field_name  = "/phiplanez"
    cmap_name   = "Blues"
    field_label = r"$\varphi/\varphi_w$"
    field_ticks = [0.0,1.0]
    xlabel      = r'$x/D$'
    
    # Limits for axes (not used in reality, but needed to use 'set_plot_settings')
    xliminf = 0.0
    xlimsup = 1.0
    yliminf = 0.0
    ylimsup = 1.0

# Streamwise vorticity
elif switcher == 1:

    Lxi         = Lz
    nxi         = nz
    field_name  = "/vortxplanex"
    cmap_name   = "RdBu"
    field_label = r"$\omega_x$"
    field_ticks = [-1.0,1.0]
    xlabel      = r'$z/D$'
    
    # Limits for axes (not used in reality, but needed to use 'set_plot_settings')
    xliminf = 0.0
    xlimsup = 1.0
    yliminf = -1.0
    ylimsup = 1.0

# Extent of the image (dimensions of the domain)
extent = [0.0, Lxi, 0.0, Ly]

#!--------------------------------------------------------------------------------------!

#!--- Mesh section and iso-levels ---!

# Create xi-coordinates vector (x or z, depending on the chosen field)
xi = np.linspace(0.0, Lxi, nxi)

# Read y-coordinates vector
y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)

#!--------------------------------------------------------------------------------------!

# Path for generic data
data_path = 'data'

# Check if the path exists and is a directory
if os.path.exists(data_path) and os.path.isdir(data_path):
   
    # Use /data to read planes if /data exists
    data_path = f"data/planes

else:

    # Use /data_ri to read planes if /data does not exists (ri specified in 'post.prm')
    data_path = f"data_r{nr:01d}/planes"
       
# Cycle on all planes
i = 0
while True:

    #!--- Reading of selected field ---!
    
    # Create the file path for the planes binary files
    file_path = data_path + field_name + f'-{i:04d}.bin' 
    
    # Exit the loop if the file does not exist
    if not os.path.exists(file_path):
        break  
    
    # Inform on current state
    print(f"We are processing file:", file_path)
    
    # Read the instantaneous field binary file into a numpy array
    with open(file_path, 'rb') as file:
        data = np.fromfile(file, dtype=np.float64)
        
    #!-------------------------------!
    
    #!--- Reading of Re_tau ---!
      
    # Create the file path for the planes .xdmf files
    file_path = data_path + field_name + f'-{i:04d}.xdmf'
    
    # Call external function to extract Re_tau value
    re_tau = extract_re_tau_value(file_path)
    
    # Print to screen the extracted value
    if re_tau is not None:
        print(f"Re_tau = {re_tau}")
    else:
        print("Re_tau value could not be extracted.")

    # Convert to string 
    re_tau = str(re_tau)

    #!--- Plotting ---!

    # Reshape the plane binary field to 2D array using Fortran order
    data = data.reshape((nxi, ny), order='F')

    # Transpose data 
    data = data.T

    # Subplots environment
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)

    # Functions to locate the colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    
    # Imshow function (unexpectedly it adjusts well the aspect ratio of the plotted image with contourf)
    im = ax.imshow(data, cmap=cmap_name, extent=extent, origin='upper')
    
    # Values of iso-levels
    lvls = np.linspace(np.min(data), np.max(data), pp.nlvl)
            
    # Plotting with filled contours    
    C = ax.contourf(xi, y, data, lvls, cmap=cmap_name)
    
    # Colorbar
    cbar = fig.colorbar(C, cax=cax, orientation='vertical', ticks=field_ticks)
       
    # Colorbar ticks 
    cbar.ax.tick_params(axis='y', labelsize=pp.fla2, length=pp.lmajt, width=pp.tick_width) 
     
    # Colorbar label
    cbar.set_label(field_label, fontsize=pp.fla, labelpad=pp.pad_cbar_lab)  

    # Axes labels and title
    ax.set_xlabel(xlabel,   fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$y/D$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_title(fr'$Re_\tau = {re_tau}$', fontsize=pp.fla)

    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

    # Specify manually the interval for xticks
    major_ticks_interval = 10.0
    ax.xaxis.set_major_locator(MultipleLocator(major_ticks_interval))

    # Saving the figure
    plt.savefig('images' + field_name + f'_{add_string}_{i:04d}.png', format='png', bbox_inches='tight', dpi=600)
    
    # Clear and close the figure to release memory
    plt.clf()
    plt.close()
    
    # Move to the next file index
    i += 1  

print(f"Reached file index {i:04d} which does not exist. Stopping.")

#!--------------------------------------------------------------------------------------!

# Change folder
#os.system("cd images")

# Create the video
#os.system("ffmpeg -y -framerate 30 -start_number 0 -c:v libx264 -i phiplanez_ttbl_small_%04d.png ttbl_phiplanez.mp4")

# Slow-down it if necessary
#os.system('fmpeg -i ttbl_phiplanez.mp4 -filter:v "setpts=4.0*PTS" ttbl_phiplanez_2.mp4')

# Delete the original
#os.system("rm ttbl_phiplanez.mp4")




