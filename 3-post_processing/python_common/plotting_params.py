#!-------------------------------------------------------------!
#! This is a small module to store Python plotting parameters. !
#!-------------------------------------------------------------!

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Settings
np.seterr(divide='ignore', invalid='ignore')

plt.rcParams.update({ 
    "text.usetex": True,  
    "font.family": "serif",
    "font.sans-serif": "Computer Modern",
    "figure.autolayout": True,
})

# Parameters for plotting
lw           = 0.6             # linewidth for plots
markersize   = 8.0             # marker size for scatter plot
fla          = 10              # fontsize of labels of x and y axes (major labels, variables) (f: fontsize, l: label, a: axes)
fla2         = 4.5             # fontsize of numbers of x and y axes 
pad_axes_lab = 2               # padding of axes labels
pad_numbers  = 3               # padding of numbers on both axes
lmajt        = 4               # length of major ticks
lmint        = 2               # length of minor ticks
tick_width   = 0.5             # width of ticks and external box

# Page settings (A4 paper format: 8.3 x 11.7 inches)
xinches      = 2.6             # size in inches in x direction of the image
yinches      = 2.2             # size in inches in y direction of the image

# Axes (line) width
mpl.rcParams['axes.linewidth'] = tick_width

# Set some useful colors
grey = [0.5, 0.5, 0.5]

# Parameter to switch between Lee & Moser reference or Cimarelli, 'Turbulence' lecture notes for the log law
iswitch = 1 # (0: Lee & Moser, 1: Cimarelli)

# Column width for writing to .txt file
c_w = 16 

# Format for numbers when written as strings (f: format, s: string)
fs = f"<{c_w}.3f"

# Format for cf only
fs2 = f"<{c_w}.8f"
