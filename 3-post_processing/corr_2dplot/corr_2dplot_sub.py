
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Subroutine used by 'corr_2plot.py'.            
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
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

# Import function to setting up, save and show plots 
from plot_subs import set_plot_settings, save_and_show_plot

#!--------------------------------------------------------------------------------------!
  
# General subroutine to plot correlation functions
def corr_2dplot(var,field_name,field_label,Lz,nz,mg_x,nu,y,ny,cmap_name,pad_cbar_lab,size_cbar,add_string,snap_numb,re_tau):
    
    # Halve the domain size and the number of points in the periodic direction to avoid periodicity effects
    # Restriction is imposed also in y if we are dealing with a Channel. 
    Lz       = Lz / 2.0                       
    nz       = nz // 2 
    var      = var[:ny,:nz]
    
    # Find the index of the maximum value
    max_index_flat = np.argmax(var)

    # Convert the flat index to a 2D index
    max_index_2d = np.unravel_index(max_index_flat, var.shape)

    """
    print(f"Index of maximum value (flat): {max_index_flat}")
    print(f"Index of maximum value (2D): {max_index_2d}")                       
    
    """
  
    # Calculate friction quantities and adimensionalize
    sh_vel   = np.sqrt(nu * np.abs(mg_x[0]))  # shear velocity (based on streamwise mean gradient)  
    delta_nu = nu / sh_vel                    # viscous length 
    y_plus   = y / delta_nu                   
    Lz_plus  = Lz / delta_nu

    # Calculate the y+ of maximum correlation
    ystar = y_plus[max_index_2d[0]] 
                                                                                   
    # Create the separation variable array (in viscous units)
    rz = np.linspace(0.0, Lz_plus, nz)

    #!--- Plot 1D section ---!

    # Limits for axes (used in 'set_plot_settings')
    xliminf = 1.0
    xlimsup = Lz_plus
    yliminf = 1.0
    ylimsup = y_plus[-1]

    # Extent of the image (dimensions of the domain)
    extent = [xliminf, xlimsup, yliminf, ylimsup]

    #!--------------------------------------------------------------------------------------!

    #!--- Mesh section and iso-levels ---!

    # Values of iso-levels        
    lvls = np.linspace(np.min(var), np.max(var), 20)

    # Calculate min and max of the array
    min_val = np.min(var)
    max_val = np.max(var)
    
    # Generate colorbar ticks as floats
    field_ticks = [min_val, max_val]

    # Format the tick labels to two decimal places
    field_tick_labels = [f"{min_val:.2f}", f"{max_val:.2f}"]
    
    # Mesh grid for contourf
    X, Y = np.meshgrid(rz, y_plus)

    #!--------------------------------------------------------------------------------------!

    # Subplots environment
    fig, ax = plt.subplots(1, 1, figsize=(pp.xinches,pp.yinches), linewidth=pp.tick_width, dpi=300)
    
    # Set the plot parameters using the function 'set_plot_settings'
    # Last argument is the switcher for semilog plot (1: yes, 0: no)
    set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, 1)
    
    # Logarithmic y-axis
    ax.set_yscale('log')

    # Functions to locate the colorbar
    divider = make_axes_locatable(ax)
    
    cax = divider.append_axes('top', size=size_cbar, pad=-0.7)
        
    # Imshow function (unexpectedly it adjusts well the aspect ratio of the plotted image with contourf)
    im = ax.imshow(var, cmap=cmap_name, extent=extent, origin='upper')
                
    # Plotting with filled contours    
    C = ax.contourf(X, Y, var, lvls, cmap=cmap_name, extend='neither')
    
    # Horizontal line to show the location of maximum correlation
    ax.hlines(y=ystar, xmin=xliminf, xmax=xlimsup, linewidth=pp.lw, color=pp.grey, linestyles='dashed')
        
    # Round y* value
    ystar = round(ystar, 1)
    
    # Text to show value of maximum correlation  
    ax.text(xlimsup*0.2, ystar*1.2, fr'$y^+ = {ystar}$', color='k', fontsize=4, ha='center')
    
    # Colorbar
    cbar = fig.colorbar(C, cax=cax, orientation='horizontal', ticks=field_ticks)
    
    """ Options commented to avoid colorbar ticks
    
    # Set colobar tick labels
    #cbar.set_ticklabels(field_tick_labels)
       
    # Colorbar ticks 
    #cbar.ax.tick_params(axis='y', labelsize=pp.fla2, length=pp.lmajt, width=pp.tick_width) 
     
    """
    
    # No ticks for the colorbar
    cbar.set_ticks([])
     
    # Colorbar label (use pp.pad_cbar_lab to use the default value for padding of the cbar label)
    cbar.set_label(field_label, fontsize=pp.fla, labelpad=pad_cbar_lab)  
    
    # Axes labels 
    ax.set_xlabel(r'$r_z^+$', fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    ax.set_ylabel(r'$y^+$',   fontsize=pp.fla, labelpad=pp.pad_axes_lab)
    
    # Save and show the figure
    save_and_show_plot(field_name, snap_numb=snap_numb, add_string=add_string, re_tau=re_tau)
    
    
    
    
