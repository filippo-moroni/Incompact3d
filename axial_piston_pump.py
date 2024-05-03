#!--------------------------------------------------------------------!
#! This script solves the continuity equation for a 0D approach       !
#! for a swash plate single-piston pump.                              !
#!                                                                    !
#! For the solution of the non-linear differential equation,          !
#! we employ backward Euler and Newton-Raphson methods.               ! 
#!--------------------------------------------------------------------!

# Libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Computer Modern Sans serif",
})

# Parameters
pi    = np.pi    # greek pi
ntot  = 3600     # total number of discretization points for the angular interval
temp  = 0.0      # temporary variable for theta calculation

cd    = 0.7      # discharge coefficient (same for delivery and suction)
omega = 3000.0   # angular velocity of the shaft of the pump [rpm]
dp    = 13.0     # piston diameter [mm]
R     = 23.0     # cylinder block diameter [mm]
beta  = 15.0     # swashplate angle [°]
V0    = 252.796  # dead volume [mm^3]
A_max = 90.0     # maximum area of suction and delivery ports [mm^2]
rho   = 850.0    # density [kg/m^3]
pd    = 100.0    # delivery pressure [bar]
ps    = 0.0      # suction  pressure [bar]
p0    = 100.0    # initial  pressure [bar]
eps   = 0.01     # tolerance for the Newton-Rapshon method
B     = 17000.0  # bulk modulus [bar]

# Conversion of units
omega = 2.0*pi*omega/60.0  # angular velocity of the shaft of the pump [rad/s]
dp    = dp / 1000.0        # piston diameter [m]
R     = 23.0 /1000.0       # cylinder block diameter [m]
beta  = np.radians(beta)   # swashplate angle [rad]
V0    = V0*(10**-9)        # dead volume [m^3]
A_max = A_max*(10**-6)     # maximum area of suction and delivery ports [m^2]
pd    = pd*(10**5)         # delivery pressure [Pa]
ps    = ps*(10**5)         # suction  pressure [Pa]
p0    = p0*(10**5)         # initial  pressure [Pa]
B     = B*(10**5)          # bulk modulus [Pa]

# Variables
theta  = np.linspace(0.0,2.0*pi,ntot)  # theta angle of rotation in radiants
dtheta = 2.0*pi/ntot                   # delta theta (angular step)

asuct   = np.zeros(ntot)               # suction area
adeliv  = np.zeros(ntot)               # delivery area
volum   = np.zeros(ntot)               # volume
volumd  = np.zeros(ntot)               # volume derivative

pp      = np.zeros(ntot)               # pressure

# Newton-Raphson method
pjp1    = 0.0                          # pressure at the j+1 iteration
pj      = 0.0                          # pressure at the j   iteration 
num     = 0.0                          # numerator of the method (g(pn+1))
den     = 0.0                          # denominator of the method (g'(pn+1))

# Backward Euler method
f       = 0.0                          # RHS of the differential equation (or for N-R method) 

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
lw = 3  # linewidth for plots
plt.plot(theta, adeliv, label="delivery", linewidth=lw) 
plt.plot(theta, asuct, label="suction", linewidth=lw)
plt.title("Delivery and suction areas", fontsize=30)
plt.xlabel(r'$\theta \,[rad]$', fontsize=30)
plt.ylabel(r'$A_s,A_d \,[m^2]$', fontsize=30)
labels = ["$0$", r"$\frac{\pi}{6}$", r"$\frac{5}{6}\pi$", r"$\pi$", r"$\frac{7}{6}\pi$", r"$\frac{11}{6}\pi$","$2 \pi$"]

# Color 'k' is black
plt.xticks([0,pi/6.0, 5.0/6.0*pi, pi, 7.0/6.0*pi, 11.0/6.0*pi, 2.0*pi], labels, color="k", size=20, rotation='horizontal')
plt.yticks(fontsize=14)
plt.tick_params(which='both', width=1)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4)

plt.legend(loc="upper left", fontsize=16)
plt.grid(which='both', color='0.65', linestyle='--', linewidth=1)
plt.show()

# Volume and volume derivative
for j in range(0,ntot):
    temp = dtheta*j
    volum[j]  = V0 + pi/4.0*(dp**2)*np.tan(beta)*R*(1-0 + np.cos(temp))
    volumd[j] = - pi/4.0*(dp**2)*np.tan(beta)*np.sin(temp)
    
# Plot volume and volume derivative
labels = ["$0$", r"$\pi$", "$2 \pi$"]

fig, ax1 = plt.subplots()
plt.title("Volume and volume derivative", fontsize=30)

color = 'tab:orange'
ax1.set_xlabel(r'$\theta \,[rad]$', fontsize=30)
ax1.set_ylabel(r'$V \,[m^3]$', fontsize=30, color=color, rotation='horizontal', ha='right')
ax1.plot(theta, volum, color=color, label="volume", linewidth=lw)
ax1.set_xticks([0, pi, 2.0*pi], labels=labels, size=20)
ax1.tick_params(axis='x', labelcolor="k", labelsize=14, which='major')
ax1.tick_params(axis='y', labelcolor=color, labelsize=14)
ax1.xaxis.grid(which='both', color='0.65', linestyle='--', linewidth=1)
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(which='both', width=1)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4)
ax1.set_xticklabels(labels=labels, fontsize=20)


# Instantiate a second axes that shares the same x-axis
ax2 = ax1.twinx()  

color = 'tab:blue'
ax2.set_ylabel(r'$\frac{d V}{d \theta} \,[m^3/rad]$', fontsize=30, color=color, rotation='horizontal', ha='left')  
ax2.plot(theta, volumd, color=color, label="volume derivative", linewidth=lw)
ax2.tick_params(axis='y', labelcolor=color, labelsize=14)

fig.tight_layout() 

lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
plt.legend(lines, labels, loc="lower right", fontsize=16)
plt.show()

## Initialize pressure
pjp1 = pd*1.5

## Backward Euler cycle
for j in range(0,ntot-1):

    temp = 100
    
    # Newton-Raphson cycle
    while temp > eps:
        
        f = B/omega/volum[j+1] * (cd * asuct [j+1] * np.sqrt(2.0/rho*np.abs(ps - pjp1)) * np.sign(ps - pjp1) + \
                                  cd * adeliv[j+1] * np.sqrt(2.0/rho*np.abs(pd - pjp1)) * np.sign(pd - pjp1) + \
                                  - omega * volumd[j+1])
    
        # Numerator of N-R  
        num = pjp1 - pp[j] - f*dtheta
        
        f = B/omega/volum[j+1] * (- cd * asuct [j+1] * 1.0 / rho / np.sqrt(2.0 / rho * np.abs(ps - pjp1)) + \
                                  - cd * adeliv[j+1] * 1.0 / rho / np.sqrt(2.0 / rho * np.abs(pd - pjp1)))
        
        # Denominator of N-R
        den = 1.0 - f*dtheta
        
        # Save this iteration of Newton-Raphson
        pj = pjp1
        
        # Newton-Raphson method
        pjp1 = pjp1 - num/den
        
        # Residual calculation
        temp = pjp1 - pj
    
    # Proceed with time integration
    f = B/omega/volum[j+1] * (cd * asuct [j+1] * np.sqrt(2.0/rho*np.abs(ps - pjp1)) * np.sign(ps - pjp1) + \
                              cd * adeliv[j+1] * np.sqrt(2.0/rho*np.abs(pd - pjp1)) * np.sign(pd - pjp1) + \
                              - omega * volumd[j+1])
    
    # Time-integration
    pp[j+1] = pp[j] + f*dtheta


## Plot of pressure
lw = 3  # linewidth for plots
plt.scatter(theta, pp, label="pressure", linewidth=lw) 
plt.title("Pressure", fontsize=30)
plt.xlabel(r'$\theta \,[rad]$', fontsize=30)
plt.ylabel(r'$P \,[Pa]$', fontsize=30)
labels = ["$0$", r"$\frac{\pi}{6}$", r"$\frac{5}{6}\pi$", r"$\pi$", r"$\frac{7}{6}\pi$", r"$\frac{11}{6}\pi$","$2 \pi$"]

# Color 'k' is black
plt.xticks([0,pi/6.0, 5.0/6.0*pi, pi, 7.0/6.0*pi, 11.0/6.0*pi, 2.0*pi], labels, color="k", size=20, rotation='horizontal')
plt.yticks(fontsize=14)
plt.tick_params(which='both', width=1)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4)

plt.legend(loc="upper left", fontsize=16)
plt.grid(which='both', color='0.65', linestyle='--', linewidth=1)
plt.show()








