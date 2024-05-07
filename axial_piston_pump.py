#!--------------------------------------------------------------------!
#! This script solves the continuity equation for a 0D approach       !
#! for a swash plate single-piston pump.                              !
#!                                                                    !
#! For the solution of the non-linear differential equation,          !
#! we employ backward Euler and Newton-Raphson.                       ! 
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
xshift_del = 0.0                       # shift for delivery port [rad] or [°]
xshift_suc = 0.0                       # shift for suction  port [rad] or [°]
q          = 0.0                       # coefficient for the straight line for ports (known term, y = mx + q)

# Conversion to [rad] if necessary
#xshift_del = np.radians(xshift_del)
#xshift_suc = np.radians(xshift_suc)

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
#fig, ax = plt.subplots()
fig, ax = plt.subplots(1, 1, figsize=(14,10))

ax.plot(theta, adeliv, label="delivery", linewidth=lw) 
ax.plot(theta, asuct, label="suction", linewidth=lw)
#plt.title("Delivery and suction areas", fontsize=30)
ax.set_xlabel(r'$\theta \,[deg]$', fontsize=50, labelpad=20)
ax.set_ylabel(r'$A_s,A_d \,[m^2]$', fontsize=50, labelpad=20)

# Limits for axes and labels
#labels = ["$0$", r"$\frac{\pi}{6}$", r"$\frac{5}{6}\pi$", r"$\pi$", r"$\frac{7}{6}\pi$", r"$\frac{11}{6}\pi$","$2 \pi$"]
labels = [r"$0^\circ$", r"$30^\circ$", r"$155^\circ$", r"$175^\circ$", r"$180^\circ$", r"$185^\circ$", r"$205^\circ$","$330^\circ$",r"$360^\circ$"]
labels = [r"$0^\circ$", r"$30^\circ$", r"$150^\circ$", r"$180^\circ$", r"$210^\circ$","$330^\circ$",r"$360^\circ$"]

# Values for ticks
v1 = pi*5.0/6.0 + xshift_del
v2 = pi + xshift_suc
v3 = pi + xshift_del
v4 = 7.0/6.0*pi + xshift_suc

#values = [0, pi/6.0, v1, v2, pi, v3, v4, 11.0/6.0*pi, 2.0*pi]
values = [0, pi/6.0, v1, pi, v4, 11.0/6.0*pi, 2.0*pi]

# Color 'k' is black
ax.set_xticks(values, labels, color="k", size=20, rotation='horizontal')
ax.set_xticklabels(labels, fontsize=24, rotation=0, ha='center')  

# Shifting to the bottom pi label to avoid superposed labels
#i = - 1
#for tick in ax.xaxis.get_major_ticks():
#    i = i + 1
#    if i == 4:
#        tick.set_pad(20)
    
ax.tick_params(axis='y', labelcolor="k", labelsize=24)

ax.tick_params(which='both', width=1)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4)

ax.legend(loc="upper left", fontsize=16)
ax.grid(which='both', color='0.65', linestyle='--', linewidth=1)

plt.savefig('delivery_suction_areas.pdf', format='pdf', bbox_inches='tight')

plt.show()

# Volume and volume derivative
for j in range(0,ntot):
    temp = dtheta*j
    volum[j]  = V0 + pi/4.0*(dp**2)*np.tan(beta)*R*(1.0 + np.cos(temp))
    volumd[j] = - pi/4.0*(dp**2)*np.tan(beta)*np.sin(temp)*R
    
# Plot volume and volume derivative
labels2 = ["$0$", r"$\pi$", "$2 \pi$"]

fig, ax1 = plt.subplots()
#plt.title("Volume and volume derivative", fontsize=30)

color = 'tab:orange'
ax1.set_xlabel(r'$\theta \,[rad]$', fontsize=30)
ax1.set_ylabel(r'$V \,[m^3]$', fontsize=30, color=color, rotation='horizontal', ha='right')
ax1.plot(theta, volum, color=color, label="volume", linewidth=lw)
ax1.set_xticks([0, pi, 2.0*pi], labels=labels2, size=20)
ax1.tick_params(axis='x', labelcolor="k", labelsize=14, which='major')
ax1.tick_params(axis='y', labelcolor=color, labelsize=14)
ax1.xaxis.grid(which='both', color='0.65', linestyle='--', linewidth=1)
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(which='both', width=1)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4)
ax1.set_xticklabels(labels=labels2, fontsize=20)


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
lines, labels2 = [sum(lol, []) for lol in zip(*lines_labels)]
plt.legend(lines, labels2, loc="lower right", fontsize=16)
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
        
    ## Time-integration, backward Euler
    pp[j+1] = pp[j] + f*dtheta    

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
ax.set_xlabel(r'$\theta \,[deg]$', fontsize=30)
ax.set_ylabel(r'$P \,[bar]$', fontsize=30)

# Color 'k' is black
ax.set_xticks(values, labels, color="k", size=20, rotation='horizontal')
ax.set_xticklabels(labels, fontsize=14, rotation=0, ha='center')  

# Shifting to the bottom pi label to avoid superposed labels
i = - 1
for tick in ax.xaxis.get_major_ticks():
    i = i + 1
    if i == 4:
        tick.set_pad(20)

# New values for y axis
#values = [min_pp, 0.0, 50.0, 100.0, max_pp]
values = [min_pp, 0.0, 50.0, 100.0]
#values = [0.0, 50.0, 100.0]

ax.set_yticks(values, size=20)
ax.tick_params(axis='y', labelcolor='k', labelsize=14, rotation=30)

# Shifting to the bottom label to avoid superposition
i = - 1
for tick in ax.yaxis.get_major_ticks():
    i = i + 1
    if i != 0:
        tick.set_pad(20)

ax.tick_params(which='both', width=1)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4)

ax.plot(min_theta,min_pp,marker='o',ms=25,mfc=(1.0,0.0,0.0,0.5),mec='None')

fig.tight_layout()
plt.legend(loc="upper right", fontsize=16)
plt.grid(which='both', color='0.65', linestyle='--', linewidth=1)
plt.show()









