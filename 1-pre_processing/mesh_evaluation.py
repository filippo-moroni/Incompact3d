"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: In this script we calculate mesh, numerics and 
!              flow parameters to setup the following flowcases:
!              - Channel flow 
!              - Temporal Turbulent Boundary Layer (TTBL)
!              See the related .pdf file for more details.
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
import sys
import os
import numpy as np

# Get the current directory
current_dir = os.path.dirname(__file__)

# Add the path to the 'python_common' directory relative to the current directory 
config_path = os.path.abspath(os.path.join(current_dir, '../', '4-common/python_common'))
sys.path.append(config_path)

# Import functions to read 'input.i3d', 'post.prm' files
from read_files import read_input_files

# Import function to setup flow parameters 
from set_flow_parameters import set_flow_parameters

# Import function to calculate stretching of the mesh
from mesh_subs import stretching_mesh_y, calculate_geometric_quantities

# Import functions to calculate memory and CPUh of simulation and to show 
# initial velocity profile
from pre_processing_tools import mem_and_cpuh, plot_initial_vel_profile

# Import fuction to create tables, add a title and write them to .txt files
from write_txt_tables import write_txt_tables

#!--------------------------------------------------------------------------------------!

# Read useful flow parameters from 'input.i3d' and 'post.prm' files
(itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, itimescheme, ioutput, ioutput_cf, iswitch_wo,  
 add_string, file1, filen, icrfile, nr, post_mean, post_vort, post_diss, post_corz, post_tke_eq
) = read_input_files('input.i3d','post.prm')

#!--- Parameters and mesh ---!
(uwall, nu, twd, y) = set_flow_parameters(re)

#!--------------------------------------------------------------------------------------!

#!--- Distinguish between flow cases (TTBL and Channel) ---!

# TTBL
if itype == 13:
    
    # Reference velocity (TTBL: wall velocity)
    uref = uwall
    
    print('>>> We are employing final BL thickness delta_99 and ') 
    print('    max cf from Cimarelli et al. (2024a),            ')
    print('    for a TTBL with final Re_tau = 500.              ')
    
    # BL thickness delta_99 of the temporal TBL at Re_tau = 500 (Cimarelli et al. (2024))  
    bl_thickness = 22.1*twd
    
    # Imposing manually that the final friction Reynolds number is 500 (only for saving to .txt file)
    re_tau = 500.0  
    
    # Maximum cf estimated at peak for TTBL (Cimarelli et al. (2024))
    cf = 0.007
    
    #!--- Ask the user to insert the domain dimensions based on TTBL thickness delta at Re_tau = 500 ---!
    print()
    scalex = float(input("Enter the scale factor for domain dimension in x-dir. based on BL thickness delta_99: "))
    print()
    
    # Common value for y-direction is 3.0 (Kozul et al. (2016), Cimarelli et al. (2024)), 
    # but values around 2.0 should still be enough to avoid confiment effects of upper boundary.
    scaley = 3.0
    
    scalez = float(input("Enter the scale factor for domain dimension in z-dir. based on BL thickness delta_99: "))
    print()
    
    # Domain dimensions
    xlx = round(scalex*bl_thickness, 3)   # domain dimension in x direction
    yly = round(scaley*bl_thickness, 3)   # domain dimension in y direction
    zlz = round(scalez*bl_thickness, 3)   # domain dimension in z direction
    
    # Ask the user the number of flow realizations
    nrealiz = int(input("Enter the number of flow realizations: "))
    print()           

# Channel    
elif itype == 3:
    
    # Reference velocity (Channel: bulk velocity)
    uref = 0.667  
    
    # Friction Reynolds number for a channel, considering the centerline Reynolds number of a laminar Poiseuille flow
    re_tau = 0.116*re**0.88
    
    print('>>> We are employing cf value from Quadrio & Ricco (2004), ') 
    print('    for a channel at Re_tau = 200.                         ')
    
    # Steady state cf of a channel flow at Re_tau = 200 (Quadrio & Ricco (2004))
    cf = 0.00793
        
    # Domain dimensions (as an example, xlx is the quantity used for calculations, while Lx is the read variable from 'input.i3d')
    xlx = Lx  # domain dimension in x direction
    yly = Ly  # domain dimension in y direction
    zlz = Lz  # domain dimension in z direction
    
    # Quadrio & Ricco (2004)
    # xlx, yly, zlz = 21.0, 2.0, 4.2
    
    # Revert to total number of ny points for channel (they are halved in 'read_input_files')
    ny = (ny - 1) * 2 + 1          
            
    # For a Channel it is supposed that we have only 1 flow realization
    nrealiz = 1

#!--- Calculations valid for both TTBL and Channel ---!

# Calculate spacings along x and z (uniform)
delta_x = xlx / nx
delta_z = zlz / nz

# Call external subroutine for stretching of the mesh
yp = stretching_mesh_y(ny, yly, beta, istret)

# Call external subroutine for printing geometric quantities of the mesh
calculate_geometric_quantities(ny, yp, delta_x)

# Shear velocity at peak cf or at steady state 
sh_vel_max = np.sqrt((cf/2.0)) * uref

# Viscous length at peak cf or at steady state
delta_nu_max = nu / sh_vel_max

# Viscous time at peak cf or at steady state
t_nu_min = nu / (sh_vel_max**2)
t_nu_min = round(t_nu_min, 4)

# First element y-dimension
delta_y1 = yp[1] - yp[0]

# Rescaling the mesh spacings with viscous unit at peak cf or at steady state
delta_x_nd_max  = round(delta_x  / delta_nu_max, 2)
delta_y1_nd_max = round(delta_y1 / delta_nu_max, 2)
delta_z_nd_max  = round(delta_z  / delta_nu_max, 2)

# Non-dimensional domain dimensions at peak cf or at steady state (nd: non-dimensional)
xlx_nd_max = round(xlx / delta_nu_max, 1)
yly_nd_max = round(yly / delta_nu_max, 1)
zlz_nd_max = round(zlz / delta_nu_max, 1)

#!--- Estimation of numerics-related parameters (CFL, D, Pé, S) at IC or at steady state ---!
   
CFL = round(uref * dt / delta_x,       2)      
D =   round(nu   * dt / (delta_y1**2), 2)    
Pe =  round(uref * delta_x / nu,       2)
        
# Stability parameter (S < 1) (see Thompson et al. (1985)) 
S = round(((uref**2)*dt)/(2.0*nu), 2)

# Calculate total number of points, number of snapshots of a single flow realization, 
# total memory requirement for all fields and estimated CPUh.
(n_tot, nsnap, mem_tot, cpuh) = mem_and_cpuh(nx,ny,nz,ifirst,ilast,itimescheme,ioutput,nrealiz)

#!--------------------------------------------------!

# Calculation of the number of mesh nodes in the viscous sublayer 
# at cf peak or at steady state.
npvis = 0     # number of points viscous sublayer
height = 0.0  # cumulative height in viscous unit (y+)
    
# Rescaling y-coordinates with viscous unit at peak cf or at steady state
yp_max = yp / delta_nu_max
      
for j in range(1, ny):
    if height + yp_max[j] - yp_max[j-1] <= 5.0:
        npvis += 1 
    height += yp_max[j] - yp_max[j-1]
        
#!--------------------------------------------------!

# This part is valid for TTBLs 
if itype == 13:
              
    # Shear velocity at Re_tau = 500:
    #sh_vel_500 = 500.0 * nu / bl_thickness
    
    # Value from G. Boga data
    sh_vel_500 = 0.04468
    
    #!--------------------------------------------------------!

    # Calculate viscous length at Re_tau = 500
    delta_nu_500  = nu / sh_vel_500
    
    # Mesh sizes at Re_tau = 500
    delta_x_nd_500  = round(delta_x  / delta_nu_500, 2)
    delta_y1_nd_500 = round(delta_y1 / delta_nu_500, 2)
    delta_z_nd_500  = round(delta_z  / delta_nu_500, 2)
                
    # Non-dimensional domain dimensions at Re_tau = 500 (nd: non-dimensional)  
    xlx_nd_500 = round(xlx / delta_nu_500, 1)
    yly_nd_500 = round(yly / delta_nu_500, 1)
    zlz_nd_500 = round(zlz / delta_nu_500, 1)
       
    #!-----------------------------------------------------------------------!
    
    # Delta y+ at the BL thickness (d: delta) at Re_tau = 500
    c = 0         # integer to store the index (c: counter)
    
    for j in range(1,ny):
        if yp[j] < bl_thickness:
            c = c + 1
    
    # We consider the element just below the estimated BL edge
    delta_yd = yp[c] - yp[c-1]
       
    # BL edge
    delta_yd_nd_500 = round(delta_yd / delta_nu_500, 1)
    
    #!---------------------------------------------!

# This part is valid for Channels 
elif itype == 3:

    # Number of points in the channel half (h: half) 
    nyh = ((ny - 1) // 2) + 1

    # Delta y+ at the channel center
    delta_yc = yp[nyh] - yp[nyh-1]
    delta_yc_nd = delta_yc / delta_nu_max
    delta_yc_nd = round(delta_yc_nd, 3)

             
#!--- Printing useful information to the screen ---!
print()
print('!----- Inputs: -----!')
print()
print('Number of mesh nodes in streamwise direction,  nx = ', nx)
print('Number of mesh nodes in wall normal direction, ny = ', ny)
print('Number of mesh nodes in spanwise direction,    nz = ', nz)
print()
print('Stretching switcher, istret          = ', istret)
print('Beta parameter, beta                 = ', beta)
print('Kinematic viscosity, nu              = ', nu)
print('Time step, dt                        = ', dt)
print('Reynolds number, Re = 1/nu           = ', re)
print('Number of flow realizations, nrealiz = ', nrealiz)
print()
print('!--- Flow case specific info ---!')
print()

if itype == 3:
    print('Reference velocity U_ref is the bulk velocity, U_bulk = ' uref)
    print()
    print('Skin friction coefficient at steady state, cf = ', cf)
    print()
    print('Estimated friction Reynolds number, Re_tau ~ ', re_tau)
    
elif itype == 13:
    print('Reference velocity U_ref is the wall velocity, Uwall = ' uref)
    print()
    print('Trip wire diameter, twd (or D) = ', twd)
    print()
    print('!--- Reference data according to Cimarelli et al. (2024): ---!')
    print()
    print('Boundary layer thickness delta_99 at Re_tau = 500, delta_99 = ', bl_thickness)
    print('Skin friction coefficient at peak, cf = ', cf)
    
print()
print('!----- Outputs: -----!')
print()
print('!--- Domain dimensions: ---!')
print()

if itype == 13:
    print('Domain dimension, Lx/D = ', xlx)     
    print('Domain dimension, Ly/D = ', yly)
    print('Domain dimension, Lz/D = ', zlz)
    print()
    print('Domain dimension at Re_tau = 500, Lx/delta_99 = ', xlx/bl_thickness)     
    print('Domain dimension at Re_tau = 500, Ly/delta_99 = ', yly/bl_thickness)
    print('Domain dimension at Re_tau = 500, Lz/delta_99 = ', zlz/bl_thickness)
    print()

if itype == 3:
    print('Domain dimension, Lx/h = ', xlx)     
    print('Domain dimension, Ly/h = ', yly)
    print('Domain dimension, Lz/h = ', zlz)

print('!--- Numerics-related parameters based on reference velocity U_ref: ---!')
print()
print('Estimated CFL,x: ', CFL)
print('Estimated D,y  : ', D)
print('Estimated Pe,x : ', Pe)
print('Estimated stability parameter S,x: ', S)
print()
print('!--- Mesh sizes at peak cf or at steady state ---!')
print()
print('Mesh size x-direction: delta_x+ = ', delta_x_nd_max)
print('Mesh size y-direction at the first element near the wall: delta_y1+ = ', delta_y1_nd_max)
print('Mesh size z-direction: delta_z+ = ', delta_z_nd_max)
print()
print('!--- Time resolution ---!')
print()
print('Estimated minimum viscous time, t_nu = ', t_nu_min)
print('Time step,                        dt = ', dt)
print('Ratio minimum viscous time and dt    = ', t_nu_min / dt)
print()

if itype == 13:
    print('!--- Mesh sizes at Re_tau = 500 (Cimarelli et al. (2024)) ---!')
    print()
    print('Mesh size x-direction: delta_x+ = ', delta_x_nd_500)
    print('Mesh size y-direction at the first element near the wall: delta_y1+ = ', delta_y1_nd_500)
    print('Mesh size z-direction: delta_z+ = ', delta_z_nd_500)
    print()
    print('Mesh size (y) at the estimated BL edge @ Re_tau = 500: delta_yd+ = ', delta_yd_nd_500)
    print()
    print('!--- Number of discretization points ---!')
    print()
    print('Number of mesh nodes in the viscous sublayer at cf peak: ', npvis)
    print()
    


elif itype == 3:
    print('Mesh size y-direction at the channel center: delta_yc+ = ', delta_yc_nd)
    
print('Total number of points: n_tot = ', n_tot)
print()
print('Number of snapshots for a single flow realization: nsnap = ', nsnap)
print()
print('Total memory requirement for snapshots [GB]: mem_tot = ', mem_tot)
print()
print('Estimated CPUh: cpuh = ', cpuh)
print()

#!-------------------------------------------------!


# table in common to both TTBL and Channel

data_input_common_1 = [
                       ["beta", "nu", "Re", "U_ref", "dt", "nrealiz", "cf", "Re_tau"],
                       [ beta,   nu,   re,   uref,    dt,   nrealiz,   cf,   re_tau ],
                      ]

data_input_common_2 = [
                       ["nx/ny/nz", "Lx/Ly/Lz" ],
                       [ nx,         xlx       ],
                       [ ny,         yly       ],
                       [ nz,         zlz       ],
                      ]  


data_output_common = [
                      ["CFL,x", "D,y", "Pe,x", "S,x"],
                      [ CFL,     D,     Pe,     S   ],
                     ]
 
data_output_common_2 = [
                        ["n_tot", "nsnap", "mem_tot [GB]", "CPUh", "sh_vel",     "t_nu",   "npvis" ],
                        [ n_tot,   nsnap,   mem_tot,        cpuh,   sh_vel_max,   t_nu_min, npvis  ],                     
                       ]



# Data only for TTBLs
if itype == 13:
 
    data_input_ttbl_1 = [
                         ["Lx,Ly,Lz/delta_99" ],
                         [ xlx/bl_thickness   ],
                         [ yly/bl_thickness   ],
                         [ zlz/bl_thickness   ],
                        ]
    
    data_input_ttbl_2 = [
                         ["delta_99 @ Re_tau = 500", "sh_vel @ Re_tau = 500" ],
                         [ bl_thickness,              sh_vel_500             ],
                        ]



    
# File creation and saving for TTBL
if itype == 13: 

    # Inside your main code, where you want to create the tables
    data_arrays = [data1, data2, data3]
    titles = ["Inputs", "Numerics-related parameters", "Outputs"]

    with open("sim_settings.txt", "w") as f:
        f.write("!----- Simulation Settings -----!\n\n")
        write_txt_tables(f, data_arrays, titles)

                       
    # Save the tables as a text file 
    with open("sim_settings.txt", "w") as f:
         f.write("!----- Temporal TBL setting parameters -----!\n")
         f.write("\n")
         f.write("!----- Inputs: -----!\n")
         f.write(table1)
         f.write("\n")
         f.write(table2)
         f.write("\n")
         f.write("!--- BL thickness delta_99 (bl_t) @ Re_tau = 500 and cf peak, according to Cimarelli et al. (2024) ---!\n")
         f.write(table3)
                 
    # Create data arrays with outputs

    data3 = [
             ["/",           "peak cf",        "Re_tau = 500"   ],
             ["delta_x+",     delta_x_nd_max,   delta_x_nd_500  ],
             ["delta_y1+",    delta_y1_nd_max,  delta_y1_nd_500 ],
             ["delta_z+",     delta_z_nd_max,   delta_z_nd_500  ],
             ["delta_yd+",   "/",               delta_yd_nd_500 ],       
            ]
           
          
    


    # Save the table as a text file and final informations
    with open("sim_settings.txt", "a") as f:
         f.write("\n")
         f.write("!----- Outputs: -----!\n")
         f.write("\n")
         f.write("!--- Non-dimensional domain dimensions: ---!\n")
         f.write(table1) 
         f.write("\n")
         f.write("!--- Numerics-related parameters: ---!\n")
         f.write(table2) 
         f.write("\n")
         f.write("!--- Mesh sizes at IC, at peak cf and at Re_tau = 500: ---!\n")
         f.write(table3) 
         f.write("\n")
         f.write("!--- Miscellaneous ---!\n")
         f.write(table4) 
         f.write("\n")
         f.write("!--- Memory (storage) requirement and CPUh ---!\n")
         f.write(table5) 
         f.write("\n")
         f.write("\n")
         f.write("!--- INFO: ---!\n")
         f.write("We are employing 3 different cf values:\n")
         f.write("\n")
         f.write("1) Value obtained from the derivative at the wall due to the IC.\n")
         f.write("\n")     
         f.write("2) Peak value found in literature (e.g. 0.007, see Cimarelli et al. (2024)).\n")
         f.write("\n")
         f.write("3) Value at Re_tau = 500, again according to Cimarelli et al. (2024).\n")
         f.write("\n")
         f.write("!--- List of acronyms & variables: ---!\n")
         f.write("\n")
         f.write("nrealiz:       Number of flow realizations considered.\n")
         f.write("S:             Stability parameter, S < 1 (Thompson et al. (1985)).\n")
         f.write("npvis:         Number of points viscous sublayer at cf peak (y+ < 5).\n")
         f.write("npsl:          Number of points initial shear layer (y+ < sl_99^+_IC).\n")
         f.write("theta_sl:      Estimated  initial momentum thickness of the shear layer (approx. 54*nu/U_wall) (dimensional).\n")
         f.write("sl_99^+_IC:    Calculated initial thickness of the shear layer (y+ where Umean < 0.01 Uwall) (non-dimensional).\n")
         f.write("sh_vel_IC:     Shear velocity of the initial condition.\n")
         f.write("sh_vel_peak:   Shear velocity at peak cf, according to Cimarelli et al. (2024).\n")
         f.write("sh_vel_500:    Shear velocity at Re_tau = 500, according to Cimarelli et al. (2024).\n")
         f.write("n_tot:         Total number of grid points.\n")
         f.write("nsnap:         Number of snapshots for a single flow realization.\n")     
         f.write("mem_tot:       Memory requirement to save snapshots in double precision, assuming 5 fields (velocity, pressure, 1 scalar field).\n")
         f.write("CPUh:          Estimated total CPUh required to complete the simulation (including different flow realizations).\n")
         f.write("               Number of elements per CPU must be higher than 100'000 and a safety factor is included.\n")
         f.write("t_nu_min:      Estimated minimum viscous time unit.\n")         
         f.write("\n")
         f.write("!-------------------------------------!\n")
         f.write("\n")
         f.write("\n")
         
# File creation and saving for Channel
elif itype == 3:
    




    # Create the tables using tabulate
    table1 = tabulate(data1, headers="firstrow", tablefmt="fancy_grid")
    table2 = tabulate(data2, headers="firstrow", tablefmt="fancy_grid")
    
    # Save the tables as a text file 
    with open("sim_settings.txt", "w") as f:
         f.write("!----- Channel setting parameters -----!\n")
         f.write("\n")
         f.write("!----- Inputs: -----!\n")
         f.write(table1)
         f.write("\n")
         f.write(table2)
    
    # Create data arrays with outputs
    data1 = [
             ["Lx+/Ly+/Lz+", "dx+/dyw+/dz+/dyc+" ],
             [ xlx_nd_max,    delta_x_nd_max     ],
             [ yly_nd_max,    delta_y1_nd_max    ],
             [ zlz_nd_max,    delta_z_nd_max     ],
             ["/",            delta_yc_nd        ],
            ]
           
    data2 = [
             ["CFL,x", "D,y", "Pe,x", "S,x"],
             [ CFL,     D,     Pe,     S   ],
            ]
                       
    data3 = [
             ["sh_vel",    "npvis", "n_tot", "nsnap", "mem_tot [GB]", "CPUh", "t_nu_min" ],
             [ sh_vel_max,  npvis,   n_tot,   nsnap,   mem_tot,        cpuh,   t_nu_min  ],
            ]
            
    # Create the tables using tabulate
    table1 = tabulate(data1, headers="firstrow", tablefmt="fancy_grid")
    table2 = tabulate(data2, headers="firstrow", tablefmt="fancy_grid")
    table3 = tabulate(data3, headers="firstrow", tablefmt="fancy_grid")
    
    # Save the table as a text file and final informations
    with open("sim_settings.txt", "a") as f:
         f.write("\n")
         f.write("!----- Outputs: -----!\n")
         f.write("\n")
         f.write("!--- Non-dimensional domain dimensions and grid spacings: ---!\n")
         f.write(table1) 
         f.write("\n")
         f.write("!--- Numerics-related parameters: ---!\n")
         f.write(table2) 
         f.write("\n")
         f.write("!--- Miscellaneous ---!\n")
         f.write(table3) 
         f.write("\n")
         f.write("!-------------------------------------!\n")
         f.write("\n")
         f.write("\n")
         f.write("!--- List of acronyms & variables: ---!\n")
         f.write("\n")
         f.write("nrealiz:       Number of flow realizations considered.\n")
         f.write("S:             Stability parameter, S < 1 (Thompson et al. (1985)).\n")
         f.write("sh_vel:        Shear velocity at steady state according to Quadrio & Ricco (2004).\n")
         f.write("npvis:         Number of points viscous sublayer at cf steady state (y+ < 5).\n")
         f.write("n_tot:         Total number of grid points.\n")
         f.write("nsnap:         Number of snapshots for a single flow realization.\n")     
         f.write("mem_tot:       Memory requirement to save snapshots in double precision, assuming 5 fields (velocity, pressure, 1 scalar field).\n")
         f.write("CPUh:          Estimated total CPUh required to complete the simulation (including different flow realizations).\n")
         f.write("               Number of elements per CPU must be higher than 100'000 and a safety factor is included.\n")
         f.write("t_nu_min:      Estimated minimum viscous time unit.\n")         
         f.write("\n")
         f.write("!-------------------------------------!\n")
         f.write("\n")
         f.write("\n") 

         
         
    


    
    



