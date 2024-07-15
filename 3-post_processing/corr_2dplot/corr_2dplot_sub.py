#!------------------------------------------------------!
#! With this script, we perform plotting of a generic   !
#! correlation function in a 2d plot as function of     !
#! y+ and spanwise separation variable rz^+.            !
#!------------------------------------------------------!

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

# Import function to setting up plots 
from plot_subs import set_plot_settings

#! General subroutine to plot correlation functions
def corr_2dplot(variable,field_name,field_label,add_string):

    # Halve the domain size and the number of points in the periodic direction to avoid periodicity effects
    # Restriction is imposed also in y if we are dealing with a Channel.
    Lz       = Lz / 2.0                       
    nz       = nz // 2 
    var      = var[:ny,:nz]                       
 
    # Calculate friction quantities and adimensionalize
    sh_vel   = np.sqrt(nu * np.abs(mg_x[0]))  # shear velocity (based on streamwise mean gradient)  
    delta_nu = nu / sh_vel                    # viscous length 
    y_plus   = y / delta_nu                   
    Lz_plus  = Lz / delta_nu
    var      = var / sh_vel
                                                                                    
    # Create the separation variable array
    rz = np.linspace(0.0, Lz, nz)
    rz = rz / delta_nu

    #!--- Plot 1D section ---!

    # Limits for axes (used in 'set_plot_settings')
    xliminf = 0.0
    xlimsup = Lz_plus
    yliminf = 0.0
    ylimsup = y_plus[-1]

    # Extent of the image (dimensions of the domain)
    extent = [xliminf, xlimsup, yliminf, ylimsup]

    #!--------------------------------------------------------------------------------------!

    #!--- Mesh section and iso-levels ---!

    # Values of iso-levels        
    lvls = np.linspace(np.min(var), np.max(var), 20)

    # Max and min values to display in the colorbar
    field_ticks = [np.min(var),np.max(var)]

    X, Y = np.meshgrid(rz, y_plus)

    #!--------------------------------------------------------------------------------------!

    # Subplots environment
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
    
    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 0)

    # Functions to locate the colorbar
    divider = make_axes_locatable(ax)
    
    cax = divider.append_axes('right', size=size_cbar, pad=0.05)
        
    # Imshow function (unexpectedly it adjusts well the aspect ratio of the plotted image with contourf)
    im = ax.imshow(var, cmap=cmap_name, extent=extent, origin='upper')
                
    # Plotting with filled contours    
    C = ax.contourf(X, Y, var, lvls, cmap=cmap_name, extend='neither')
    
    # Colorbar
    cbar = fig.colorbar(C, cax=cax, orientation='vertical', ticks=field_ticks)
       
    # Colorbar ticks 
    cbar.ax.tick_params(axis='y', labelsize=pp.fla2, length=pp.lmajt, width=pp.tick_width) 
     
    # Colorbar label (use pp.pad_cbar_lab to use the default value for padding of the cbar label)
    cbar.set_label(field_label, fontsize=pp.fla, labelpad=pad_cbar_lab)  
    
    # Axes labels 
    ax.set_xlabel(r'$r_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$y^+$',   fontsize=pp.fla, labelpad=pp.pad_axes_lab)

    # Show the plot
    plt.show()
    
    # Saving the figure
    plt.savefig('plots' + field_name + f'_{add_string}.pdf', format='pdf', bbox_inches='tight', dpi=600)
    
    # Clear and close the figure to release memory
    plt.clf()
    plt.close()
    
    
