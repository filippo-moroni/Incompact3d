#!-------------------------------------------------------!
#! This is a small function to setting up semilog plots. !
#!-------------------------------------------------------!

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

def set_semilog_plot_settings(ax, xliminf, xlimsup, yliminf, ylimsup, pp):

    # Axes limits
    plt.xlim([xliminf, xlimsup])
    plt.ylim([yliminf, ylimsup])

    # Logarithmic x-axis and linear y-axis
    ax.set_xscale('log')
    ax.set_yscale('linear')

    # Minor x-ticks based on log10
    ax.xaxis.set_minor_locator(LogLocator(base=10,subs='all'))
    
    # Setting major and minor ticks on both axes
    ax.tick_params(axis='both', which='major', direction='in', length=pp.lmajt, width=pp.tick_width, pad=pp.pad_numbers, labelsize=pp.fla2, labelcolor='k') 
    ax.tick_params(axis='both', which='minor', direction='in', length=pp.lmint, width=pp.tick_width)

    # Setting x-ticks labels
    plt.xticks(ha='left')
