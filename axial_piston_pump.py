#!--------------------------------------------------------------------!
#! This script solves the continuity equation for a 0D approach       !
#! for a swash plate single-piston pump.                              !
#!                                                                    !
#! For the solution of the non-linear implicit differential equation, !
#! we employ backward Euler and Newton-Raphson methods.               ! 
#!--------------------------------------------------------------------!

# Libraries
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

# Parameters
A_max = 90.0e-6  # maximum area of suction and delivery ports [m^2]
pi    = np.pi    # greek pi
ntot  = 720      # total number of discretization points for the angular interval
temp  = 0.0      # temporary variable for theta calculation

omega = 3000     # rpm of the shaft of the pump
omega = 2.0*pi*omega/60.0  # conversion to radiants [rad/s]

cd = 0.7         # discharge coefficient (same for delivery and suction)

# Variables
theta  = np.linspace(0.0,2.0*pi,ntot)  # theta angle of rotation in radiants
dtheta = 2.0*pi/ntot                   # delta theta (angular step)

asuct   = np.zeros(ntot)               # suction area
adeliv  = np.zeros(ntot)               # delivery area

# Delivery area
for j in range(0,ntot):
    temp = dtheta*j
    if temp < pi/6.0:
        adeliv[j] = A_max/(pi/6.0) * j * dtheta
    elif temp >= pi/6.0 and temp < pi*5.0/6.0:
        adeliv[j] = A_max
    elif temp >= pi*5.0/6.0 and temp < pi:
        adeliv[j] = - A_max/(pi/6.0) * j * dtheta + 6.0*A_max
            
# Suction area
for j in range(0,ntot):
    temp = dtheta*j
    if temp >= pi and temp < 7.0/6.0*pi:
        asuct[j] = A_max/(pi/6.0) * j * dtheta - 6.0*A_max
    elif temp >= 7.0/6.0*pi and temp < 11.0/6.0*pi:
        asuct[j] = A_max
    elif temp >= 11.0/6.0*pi and temp <= 2.0*pi:
        asuct[j] = - A_max/(pi/6.0) * j * dtheta + 12.0*A_max
                     
# Plot of delivery and suction areas
plt.scatter(theta, adeliv, label="delivery") 
plt.scatter(theta, asuct, label="suction")
plt.title("Delivery and suction areas", fontsize=30)
plt.xlabel(r'$\theta$', fontsize=30)
plt.ylabel(r'$A_s,A_d$', fontsize=30)
labels = ["$0$", r"$\pi$", "$2 \pi$"]

# Color 'k' is black
plt.xticks([0, pi, 2.0*pi], labels, color="k", size=20, rotation='horizontal')
plt.legend(loc="upper left", fontsize=16)
plt.show()

# Volume













