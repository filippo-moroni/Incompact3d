
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: This is a small module to store Python plotting parameters and
!              constans (e.g. double precision numbers).
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

"""
!----------------------------------------! 
! Parameters and libraries for plotting.
!----------------------------------------!
"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings

# Set Numpy
np.seterr(divide='ignore', invalid='ignore')

# Set Latex 
plt.rcdefaults() 
plt.rcParams.update({ 
    "text.usetex": True,  
    "font.family": "serif",
    "font.serif": "Computer Modern",
    "figure.autolayout": True,
})
     
# Filter warnings
warnings.filterwarnings("error")

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

"""
# Format for numbers when written as strings (f: format, s: string)
(>: right-adjust, <: left-adjust).
"""
fs  = f">{c_w}.3f"   # 3 decimal digits
fs6 = f">{c_w}.6f"   # 6 decimal digits
fs8 = f">{c_w}.8f"   # 8 decimal digits (used mainly for friction coefficient cf)

# Parameters for contourf plots
nlvl = 256           # number of isolevels
pad_cbar_lab = -8    # padding of label of colorbar (cbar)

#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

"""
!--------------------------------------------! 
! Parameters and libraries for calculations.
!--------------------------------------------!
"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings

# Set Numpy
np.seterr(divide='ignore', invalid='ignore')



