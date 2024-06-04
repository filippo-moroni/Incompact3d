#!----------------------------------------------!
#! With this script we plot the correlation     !
#! functions calculated from 'post_incompact3d' !
#!----------------------------------------------!

import numpy as np
import matplotlib.pyplot as plt

# Open the file in binary mode
with open('data_post/corr_stats.bin', 'rb') as f:
    # Read the data into a NumPy array
    corr_z = np.fromfile(f, dtype=np.float64)       # Change dtype according to your data

#print (len(corr_z))

corr_z = np.reshape(corr_z, (64, 65, 3), order='F') # Fortran-like index ordering

Ruuz   = corr_z[:,3,0] 

rz = np.linspace(0, 4.0, 64)

plt.plot(rz, Ruuz)
plt.show()
