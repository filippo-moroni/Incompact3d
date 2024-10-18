
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this script, we perform saving of images from 
!              instantaneous planes of scalar field or of streamwise   
!              vorticity.                                                                                         
!              Inspired by 'snap2png.py' by R. Corsini.            
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
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
config_path = os.path.abspath(os.path.join(current_dir, '../../4-common', 'python_common'))
sys.path.append(config_path)

# Import the plotting_params module
import plot_params as pp

# Import function to set plots
from plot_subs import set_plot_settings

# Import function to read 'input.i3d' and 'post.prm' files
from read_files import read_input_files

# Import function to read friction Reynolds number Re_tau from .xdmf files
from ttbl_subs import extract_re_tau_value 

#!--------------------------------------------------------------------------------------!

# Print to screen what the program does

print("!--- 'field_images.py' ---!")
print()
print(" Creation and saving of .png images from instantaneous")
print(" planes of scalar field or of streamwise vorticity.")
print()
print(" It has a recursive nature, thus it can be used to create")
print(" several images for creation of videos for flow visualization.")
print()

#!--------------------------------------------------------------------------------------!

# Create folders to store images
os.makedirs('images', mode=0o777, exist_ok=True)

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
(itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, itimescheme, ioutput, ioutput_cf, iswitch_wo,  
 add_string, file1, filen, icrfile, nr, post_mean, post_grad, post_corz, post_tke_eq
) = read_input_files('input.i3d','post.prm')

#!--------------------------------------------------------------------------------------!

# Asking the user what he wants to plot (scalar field, streamwise vorticity)
print()
switcher = int(input(">>> Select a field (0: scalar field, 1: streamwise vorticity): "))
print()

# Asking the user the realization folder to use (if TTBL)
if itype == 13:

    realiz = int(input(">>> Specify the realization folder to use: "))
    print()

# Asking the user if he wants to plot the title (that is friction Reynolds number)
i_title = int(input(">>> Add title to .png images? (Friction Reynolds number) (0: no, 1: yes): "))
print()

# Scalar field
if switcher == 0:

    Lxi           = Lx
    nxi           = nx
    field_name    = "/phiplanez"
    cmap_name     = "Blues"
    field_label   = r"$\varphi/\varphi_w$"
    field_ticks   = [0.0,1.0]
    xlabel        = r'$x/D$'
    pad_cbar_lab  = -8
    size_cbar     = '5%'
    maj_ticks_int = 10.0
        
# Streamwise vorticity
elif switcher == 1:

    Lxi           = Lz
    nxi           = nz
    field_name    = "/vortxplanex"
    cmap_name     = "RdBu"
    field_label   = r"$\omega_x D/U_w$"
    field_ticks   = [-0.5,0.5]
    xlabel        = r'$z/D$'
    pad_cbar_lab  = -18
    size_cbar     = '5%'
    maj_ticks_int = 10.0
    
# Create subfolder of the selected field 
os.makedirs('images' + field_name, mode=0o777, exist_ok=True)
    
# Extent of the image (dimensions of the domain)
# (if too large, we can try to rescale through a factor that multiplies the array)
extent = [0.0, Lxi, 0.0, Ly]

# Limits for axes (used in 'set_plot_settings')
xliminf = 0.0
xlimsup = Lxi
yliminf = 0.0
ylimsup = Ly

#!--------------------------------------------------------------------------------------!

#!--- Mesh section and iso-levels ---!

# Create xi-coordinates vector (x or z, depending on the chosen field)
xi = np.linspace(0.0, Lxi, nxi)

# Read y-coordinates vector
y = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)

# Values of iso-levels        
lvls = np.linspace(field_ticks[0], field_ticks[1], pp.nlvl)

#!--------------------------------------------------------------------------------------!

# Path for generic data
data_path = 'data'

# Check if the path exists and is a directory
if os.path.exists(data_path) and os.path.isdir(data_path):
   
    # Use /data to read planes if /data exists
    data_path = f"data/planes"

else:

    # Use /data_r1 to read planes if /data does not exists 
    data_path = f"data_r{realiz}/planes"
       
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
    print(f">>> We are processing file:", file_path)
    
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
        print(f">>> Re_tau = {re_tau}")
    else:
        print(">>> Re_tau value could not be extracted.")

    # Convert to string 
    re_tau = str(re_tau)

    #!--- Plotting ---!

    # Reshape the plane binary field to 2D array using Fortran order
    
    # Scalar field
    if switcher == 0:
        data = data.reshape((nxi, ny), order='F')
        
        # Transpose data 
        data = data.T
              
    # Streamwise vorticity
    elif switcher == 1:
        data = data.reshape((ny, nxi), order='F')
            
    # Modify the array of data to limit to the specified range
    data = np.where((data < field_ticks[0]), field_ticks[0], data)
    data = np.where((data > field_ticks[1]), field_ticks[1], data)
                   
    # Subplots environment
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
    
    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

    # Functions to locate the colorbar
    divider = make_axes_locatable(ax)
    
    cax = divider.append_axes('right', size=size_cbar, pad=0.05)
        
    # Imshow function (unexpectedly it adjusts well the aspect ratio of the plotted image with contourf)
    im = ax.imshow(data, cmap=cmap_name, extent=extent, origin='upper')
                
    # Plotting with filled contours    
    C = ax.contourf(xi, y, data, lvls, cmap=cmap_name, extend='neither')
    
    # Colorbar
    cbar = fig.colorbar(C, cax=cax, orientation='vertical', ticks=field_ticks)
       
    # Colorbar ticks 
    cbar.ax.tick_params(axis='y', labelsize=pp.fla2, length=pp.lmajt, width=pp.tick_width) 
     
    # Colorbar label (use pp.pad_cbar_lab to use the default value for padding of the cbar label)
    cbar.set_label(field_label, fontsize=pp.fla, labelpad=pad_cbar_lab)  
    
    # Axes labels
    ax.set_xlabel(xlabel,   fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$y/D$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    # Title
    if i_title == 1:
        ax.set_title(fr'$Re_\tau = {re_tau}$', fontsize=pp.fla)

    # Specify manually xticks
    ax.xaxis.set_major_locator(MultipleLocator(maj_ticks_int))

    # Saving the figure
    plt.savefig('images' + field_name + field_name + f'_{add_string}_{i:04d}.png', format='png', bbox_inches='tight', dpi=600)
    
    # Clear and close the figure to release memory
    plt.clf()
    plt.close()
    
    # Move to the next file index
    i += 1  

print(f">>> Reached file index {i:04d} which does not exist. Stopping.")

#!--------------------------------------------------------------------------------------!

# Change folder
#os.system("cd images")

# Create the video
#os.system("ffmpeg -y -framerate 30 -start_number 0 -i field_flowcasename_%04d.png field_flowcasename.mp4")

#!--------------------------------------------------------------------------------------!




