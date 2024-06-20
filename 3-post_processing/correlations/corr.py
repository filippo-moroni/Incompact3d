#!----------------------------------------------!
#! With this script we plot the correlation     !
#! functions calculated from 'post_incompact3d' !
#!----------------------------------------------!

import numpy as np
import matplotlib.pyplot as plt

# Read input file
with open('input.i3d', 'r') as file:
    
    # Read all lines into a list
    lines = file.readlines()
    
    # Extract itype, nx, nz, Lx, Ly, Lz, Re 
    itype = lines[7]  
    nx    = lines[14]
    nz    = lines[16]
    Lx    = lines[22]
    Ly    = lines[23]
    Lz    = lines[24]
    re    = lines[27]
    
    # Removing characters in front of the extracted strings and the comments
    itype = itype.split('=')[-1].strip()
    
    nx    = nx.split('!')[0]
    nx    = nx.split('=')[-1].strip()
    
    nz    = nz.split('!')[0]
    nz    = nz.split('=')[-1].strip()
    
    Lx    = Lx.split('!')[0]
    Lx    = Lx.split('=')[-1].strip()
    
    Ly    = Ly.split('!')[0]
    Ly    = Ly.split('=')[-1].strip()
    
    Lz    = Lz.split('!')[0]
    Lz    = Lz.split('=')[-1].strip()
    
    re    = re.split('!')[0]
    re    = re.split('=')[-1].strip()
    
    # Convert to integer
    itype = int(itype)
    nx    = int(nx)
    nz    = int(nz)
    Lx    = np.float64(Lx)
    Ly    = np.float64(Ly)
    Lz    = np.float64(Lz)
    re    = np.float64(re)

#!--- Reading of files section ---!
print()

# Channel
if itype == 3:
    
    print("!--- Plotting of correlations for a channel ---!")

    # Reading of mean statistics
    M1 = np.loadtxt('data_post/corr_stats.txt', skiprows=0, delimiter=None, dtype=np.float64)
    
# TTBL
elif itype == 13:

    print("!--- Plotting of correlations for a TTBL ---!")

    # Asking to the user the specific snapshot to show
    snap_numb = input("Enter the snapshot number to show (4 digits): ")
         
    # Reading of mean statistics
    M1 = np.loadtxt(f'data_post/corr_stats-{snap_numb}.txt', skiprows=0, delimiter=None, dtype=np.float64)
    
print()

# Reading of grid points
y = np.loadtxt('yp.dat')

# Number of points in y direction
ny = len(y)

# Halve the points in y direction for a channel
if itype == 3:
    ny = (ny - 1) // 2 + 1

# Select the height at which correlations are plotted
yplus_in = input("Enter y+ value of the plane: ")


# Extracting quantities from the full matrix
mean_u  = M1[:,0]   
var_u   = M1[:,3]


# Open the file in binary mode
with open('data_post/corr_stats.bin', 'rb') as f:
    # Read the data into a NumPy array
    corr_z = np.fromfile(f, dtype=np.float64)       # Change dtype according to your data

#print (len(corr_z))

corr_z = np.reshape(corr_z, (64, 65, 3), order='F') # Fortran-like index ordering

Ruuz   = corr_z[:,0,0] 

rz = np.linspace(0, 4.0, 64)

plt.plot(rz, Ruuz)
plt.show()
