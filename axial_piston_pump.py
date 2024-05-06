#!--------------------------------------------------------------------!
#! This script solves the continuity equation for a 0D approach       !
#! for a swash plate single-piston pump.                              !
#!                                                                    !
#! For the solution of the non-linear differential equation,          !
#! we employ BDF3 and Newton-Raphson.                                 ! 
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
ntot  = 7200     # total number of discretization points for the angular interval [0, 2*pi]
temp  = 0.0      # temporary variable 

cd    = 0.7      # discharge coefficient (same for delivery and suction)
omega = 3000.0   # angular velocity of the shaft of the pump [rpm]
dp    = 13.0     # piston diameter [mm]
R     = 23.0     # cylinder block diameter [mm]
beta  = 15.0     # swashplate angle [°]
V0    = 252.796  # dead volume [mm^3]
A_max = 90.0     # maximum area of suction and delivery ports [mm^2]
rho   = 850.0    # density of work fluid [kg/m^3]
pd    = 100.0    # delivery pressure [bar]
ps    = 0.0      # suction  pressure [bar]
p0    = 100.0    # initial  pressure [bar]
eps   = 0.001    # tolerance for the Newton-Raphson method
B     = 15000.0  # bulk modulus of work fluid [bar]
huge  = 10**9    # a generic large number 

# Conversion of units
omega = 2.0*pi*omega/60.0  # angular velocity of the shaft of the pump [rad/s]
dp    = dp / 1000.0        # piston diameter [m]
R     = R / 1000.0         # cylinder block diameter [m]
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
pjp1    = np.double(0.0)               # pressure at the j+1 iteration
pj      = np.double(0.0)               # pressure at the j   iteration 
num     = np.double(0.0)               # numerator of the method (g(pn+1))
den     = np.double(0.0)               # denominator of the method (g'(pn+1))

# Numerics 
f       = np.double(0.0)               # RHS of the differential equation or function of NR method
fprime  = np.double(0.0)               # First derivative for NR method

# Shifts for ports 
xshift_del =  -2.0                  # shift for delivery port [°]
xshift_suc =   2.0                  # shift for suction  port [°]
q          = 0.0                       # coefficient for the straight line for ports (known term, y = mx + q)

xshift_del = np.radians(xshift_del)
xshift_suc = np.radians(xshift_suc)

# Delivery area
for j in range(0,ntot):

    temp = dtheta*j
    
    if temp < pi/6.0:
        adeliv[j] = A_max/(pi/6.0) * j * dtheta
        
    elif temp >= pi/6.0 and temp < pi*5.0/6.0 + xshift_del:
        adeliv[j] = A_max
        
    elif temp >= pi*5.0/6.0 + xshift_del and temp < pi + xshift_del:
        q = 6.0*A_max/pi*(pi + xshift_del)
        adeliv[j] = - A_max/(pi/6.0) * j * dtheta + q
            
# Suction area
for j in range(0,ntot):

    temp = dtheta*j
    
    if temp >= pi + xshift_suc and temp < 7.0/6.0*pi + xshift_suc:
        q = -6.0*A_max/pi*(pi + xshift_suc)
        asuct[j] = A_max/(pi/6.0) * j * dtheta + q
        
    elif temp >= 7.0/6.0*pi + xshift_suc and temp < 11.0/6.0*pi:
        asuct[j] = A_max
        
    elif temp >= 11.0/6.0*pi and temp <= 2.0*pi:
        asuct[j] = - A_max/(pi/6.0) * j * dtheta + 12.0*A_max
                     
# Plot of delivery and suction areas
lw = 3  # linewidth for plots
plt.plot(theta, adeliv, label="delivery", linewidth=lw) 
plt.plot(theta, asuct, label="suction", linewidth=lw)
#plt.title("Delivery and suction areas", fontsize=30)
plt.xlabel(r'$\theta \,[rad]$', fontsize=30)
plt.ylabel(r'$A_s,A_d \,[m^2]$', fontsize=30)
labels = ["$0$", r"$\frac{\pi}{6}$", r"$\frac{11}{12}\pi$", r"$\pi$", r"$\frac{13}{12}\pi$", r"$\frac{11}{6}\pi$","$2 \pi$"]

# Color 'k' is black
plt.xticks([0, pi/6.0, 11.0/12.0*pi, pi, 13.0/12.0*pi, 11.0/6.0*pi, 2.0*pi], labels, color="k", size=20, rotation='horizontal')
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
    volum[j]  = V0 + pi/4.0*(dp**2)*np.tan(beta)*R*(1.0 + np.cos(temp))
    volumd[j] = - pi/4.0*(dp**2)*np.tan(beta)*np.sin(temp)*R
    
# Plot volume and volume derivative
labels = ["$0$", r"$\pi$", "$2 \pi$"]

fig, ax1 = plt.subplots()
#plt.title("Volume and volume derivative", fontsize=30)

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
ax2.tick_params(which='both', width=1)
ax2.tick_params(which='major', length=7)
ax2.tick_params(which='minor', length=4)

fig.tight_layout() 

lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
plt.legend(lines, labels, loc="lower right", fontsize=16)
plt.show()


## Initialize pressure
pp[:] = p0

## Backward Euler cycle
for j in range(0,ntot-1):

    # Initial temporary value for the residual
    temp = huge
    
    # Initial guess for Newton-Raphson method
    pjp1 = pp[j] * 2.0
    
    # Newton-Raphson cycle, with exit condition with tolerance equal to 'eps'
    while temp > eps:
        
        f = (B/omega/volum[j+1] * (cd * asuct [j+1] * np.sqrt(2.0/rho*np.abs(ps - pjp1)) * np.sign(ps - pjp1) + \
                                +  cd * adeliv[j+1] * np.sqrt(2.0/rho*np.abs(pd - pjp1)) * np.sign(pd - pjp1) + \
                                - omega * volumd[j+1]))
    
        # Numerator of N-R  
        num = pjp1 - pp[j] - f*dtheta
        

        fprime = (- B/omega/volum[j+1] * (cd * asuct [j+1] / np.sqrt(2.0 * rho * np.abs(ps - pjp1))   + \
                                       +  cd * adeliv[j+1] / np.sqrt(2.0 * rho * np.abs(pd - pjp1))))
        
        # Denominator of N-R
        den = 1.0 - fprime*dtheta
        
        # Store of current iteration
        #pj = pjp1
                
        # Newton-Raphson method
        pjp1 = pjp1 - num/den
        
        # New estimation of f
        f = (B/omega/volum[j+1] * (cd * asuct [j+1] * np.sqrt(2.0/rho*np.abs(ps - pjp1)) * np.sign(ps - pjp1) + \
                                +  cd * adeliv[j+1] * np.sqrt(2.0/rho*np.abs(pd - pjp1)) * np.sign(pd - pjp1) + \
                                - omega * volumd[j+1]))
        
        # Residual calculation (total function of N-R to be forced to zero)
        temp = pjp1 - pp[j] - f*dtheta
        
        # We can also enforce the residual on the pressure calculation itself (relative residual)
        #temp = (pjp1 - pj)/pjp1
        
        # Absolute value for the residual
        temp = np.abs(temp)
        
    ## Time-integration
    
    # Start with backward Euler
    if j == 0:
        pp[j+1] = pp[j] + f*dtheta
    
    # Backward differencing 2nd order (BDF2)    
    elif j == 1:
        pp[j+1] = 4.0/3.0*pp[j] - 1.0/3.0*pp[j-1] + 2.0/3.0*f*dtheta
    
    # Backward differencing 3rd order (BDF3)    
    elif j > 1:
        pp[j+1] = 18.0/11.0*pp[j] - 9.0/11.0*pp[j-1] + 2.0/11.0*pp[j-2] + 6.0/11.0*f*dtheta
    

# Rescaling of pressure for plotting
pp = pp/(10**5)  # [bar]

# Maximum and minimum of pressure
max_pp = max(pp)
min_pp = min(pp)

# Angle theta for the minimum pressure
index = np.where(pp == min_pp)
min_theta = theta[index]

## Plot of pressure
fig, ax = plt.subplots()

ax.scatter(theta, pp, label="pressure", marker="o", s=20, facecolors='C0', edgecolors='C0') 
#plt.title("Pressure", fontsize=30)
ax.set_xlabel(r'$\theta \,[rad]$', fontsize=30)
ax.set_ylabel(r'$P \,[bar]$', fontsize=30)

# Limits for axes and labels
labels = ["$0$", r"$\frac{\pi}{6}$", r"$\frac{11}{12}\pi$", r"$\pi$", r"$\frac{13}{12}\pi$", r"$\frac{11}{6}\pi$","$2 \pi$"]
values = [0,pi/6.0, 11.0/12.0*pi, pi, 13.0/12.0*pi, 11.0/6.0*pi, 2.0*pi]
ax.set_xticks(values, labels=labels, size=20)


values = [min_pp, 0.0, 50.0, 100.0, max_pp]
ax.set_yticks(values, size=20)

ax.tick_params(axis='y', labelcolor='k', labelsize=14)
ax.tick_params(which='both', width=1)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4)

#ax.scatter(min_theta,min_pp,marker='o',s=20,facecolors='none', edgecolors='r')
ax.plot(min_theta,min_pp,marker='o',ms=20,mfc=(1.0,0.0,0.0,0.5),mec='None')

fig.tight_layout()
plt.legend(loc="lower left", fontsize=16)
plt.grid(which='both', color='0.65', linestyle='--', linewidth=1)
plt.show()








