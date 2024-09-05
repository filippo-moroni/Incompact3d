#!-----------------------------------------------!
#! In this file, we store useful subroutines     !
#! for plotting with Python matplotlib.          !
#!-----------------------------------------------!

#!-----------------------------------------------!
#! This is a small function to setting up plots. !
#!-----------------------------------------------!

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

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

#!-------------------------------------------------!
#! With this script, we close the plot section for !
#! statistics plotting.                            !
#!-------------------------------------------------! 

import matplotlib.pyplot as plt

def save_and_show_plot(variable_name, snap_numb=None, add_string=None, re_tau=None, y_plus_in=None, subfolder=None, description=None):
    
    """
    Saves and shows a plot with a given variable name and optional parameters.
    
    Parameters:
    - variable_name (str):           The name of the variable to be saved in the filename.
    - snap_numb (str, optional):     Snapshot number to be included in the filename. Added if itype is 13 (TTBL).
    - add_string (str, optional):    Additional string to be included in the filename, used to add the flowcase name.
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
    


