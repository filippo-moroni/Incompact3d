#!-----------------------------------------------!
#! This is a small function to setting up plots. !
#!-----------------------------------------------!

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

def set_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp, iswitch_slp):

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
    

