
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: In this file, we store useful subroutines for plotting with 
!              Python matplotlib.                         
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Common libraries
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '..', 'python_common'))
sys.path.append(config_path)

# Import the plotting_params module
import plot_params as pp

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: This is a small function to setting up plots. 
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, iswitch_slp):

    """
    Parameters:
    
    iswitch_slp: Switcher to enable semilog plot for x-axis (0: no, 1: yes).
        
    """

    # Axes limits
    plt.xlim([xliminf, xlimsup])
    plt.ylim([yliminf, ylimsup])

    # Axes scales definition
    if iswitch_slp == 1:
    
        # Logarithmic x-axis
        ax.set_xscale('log')
        
        # Minor x-ticks based on log10
        ax.xaxis.set_minor_locator(LogLocator(base=10,subs='all'))
    
    elif iswitch_slp == 0:
    
        # Linear x-axis   
        ax.set_xscale('linear')
                    
    # Linear y-axis
    ax.set_yscale('linear')

    # Setting major and minor ticks on both axes
    ax.tick_params(axis='both', which='major', direction='in', length=pp.lmajt, width=pp.tick_width, pad=pp.pad_numbers, labelsize=pp.fla2, labelcolor='k') 
    ax.tick_params(axis='both', which='minor', direction='in', length=pp.lmint, width=pp.tick_width)

    # Setting x-ticks labels
    plt.xticks(ha='left')
 
#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this script, we close the plot section for statistics 
!              plotting (save and show).
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def save_and_show_plot(variable_name, snap_numb=None, ts=None, add_string=None, re_tau=None, y_plus_in=None, subfolder=None, description=None):
    
    """
    Saves and shows a plot with a given variable name and optional parameters.
    
    Parameters:
    - variable_name (str):           The name of the variable to be saved in the filename.
    - snap_numb (str, optional):     Snapshot number to be included in the filename. Added if itype is 13 (TTBL).
    - ts (str, optional):            Time-step at which we are plotting. Used only if statistics from 'cf_monitoring' are plotted.
    - add_string (str, optional):    Additional string to be included in the filename, used to add the flowcase name (add: additional).
    - y_plus_in (float64, optional): The additional parameter to be included in the filename if we are plotting and saving correlation functions.
    - re_tau (int, optional):        Additional Friction Reynolds number value to add the the filename if we are saving a plot of a TTBL.
    - subfolder (str, optional):     Name of the subfolder in which we are saving the .pdf file.
    - description (str, optional):   Description to be added to the PDF metadata.
    
    """
    
    # Initialize the filename
    filename = f'plots'
      
    # Add subfolder if it is provided
    if subfolder is not None:
        filename += f'/{subfolder}'
        
    # Add variable name
    filename += f'/{variable_name}'
    
    # Add snap_numb if it is provided
    if snap_numb is not None:
        filename += f'-{snap_numb}'
    
    # Add time-step if it is provided
    if ts is not None:
        filename += f'-ts{ts}'
    
    # Add add_string if it is provided
    if add_string is not None:
        filename += f'_{add_string}'
  
    # Add re_tau if it is provided
    if re_tau is not None:
        filename += f'_retau={re_tau}'
    
    # Add y+ location of correlations if provided
    if y_plus_in is not None:
        filename += f'_y+={y_plus_in}'
        
    # Add .pdf extension to filename
    filename += '.pdf'
    
    # Add metadata with supported fields
    metadata = {'Title': variable_name,
                'Author': 'Filippo Moroni',
                'Subject': description}  
    
    # Save the figure and plot it
    plt.savefig(filename, format='pdf', bbox_inches='tight', dpi=600, metadata=metadata)
    plt.show()

# Example usage:
# save_and_show_plot('umean', snap_numb=snap_numb, add_string=add_string, re_tau=re_tau)
# save_and_show_plot('Cuuz', add_string=add_string, y_plus_in=y_plus_in)
    
#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this script, we setup the mean streamwise velocity
!              profile as reference for plotting. Different log law constants
!              are employed depending on the flowcase.
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

def get_ref_mean_vel_profile(itype,iswitch):

    # Viscous sub-layer
    y_plus_vsl = np.linspace(1, 15, 15)
    u_plus_vsl = y_plus_vsl

    # Log law constants based on specific flow case
    if itype == 13:
   
        # Kozul et al. (2016)
        k = 0.384
        B = 4.173
        
    elif itype == 3:

        if iswitch == 0:
    
            # Lee and Moser (2015)
            k = 0.384
            B = 4.27
    
        elif iswitch == 1:
        
            # Cimarelli ('Turbulence' lecture notes)
            k = 0.37
            B = 5.2

    # Von Karman law
    y_plus_k = np.linspace(5, 180, 175)
    u_plus_k = (1.0 / k) * np.log(y_plus_k) + B

    return (y_plus_vsl,u_plus_vsl,y_plus_k,u_plus_k)
    
#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!


