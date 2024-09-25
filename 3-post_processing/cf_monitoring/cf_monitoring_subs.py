
"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this subroutine we perform calculation of 6th order 
!              accurate calculations of integrals of TTBL thickness 
!              parameters using 'umean' files printed at the same time as 
!              'cf_monitoring'.  
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
from scipy.interpolate import InterpolatedUnivariateSpline

def calculate_thickness_param():

    #!--- Reading of files section and setup of flow parameters ---!

    # Read useful flow parameters from 'input.i3d' and 'post.prm' files
    (itype, nx, ny, nz, istret, beta, Lx, Ly, Lz, re, dt, ifirst, ilast, numscalar, ioutput, ioutput_cf, iswitch_wo,  
     add_string, file1, filen, icrfile, nr, post_mean, post_vort, post_diss, post_corz, post_tke_eq
    ) = read_input_files('input.i3d','post.prm')

    # Reading of yp coordinates
    yp = np.loadtxt('yp.dat', delimiter=None, dtype=np.float64)

    # First element of yp vector (y = 0)
    y0 = yp[0]

    # Last  element of yp vector (y = Ly, height of the domain)    
    yn = yp[-1]

    #!--- Parameters ---!
    uwall, nu, twd = set_flow_parameters(itype, re)
                          
    #!--------------------------------------------------------------------------------------!

    # First time-step is skipped at the moment for the reading of umean
    ts1 = ioutput_cf
                
    # Last time-step
    tsn = ilast
    
    # Number of savings due to 'print_cf' subroutine of Incompact3d solver 'modified'
    nsavings = (tsn - ts1) // ioutput_cf + 1

    # Initialize the mean streamwise velocity profile array (function of y and specific flow realization)
    umean = np.zeros(ny, nr)

    # Initialize the array to sum the mean streamwise velocity profile (function of y and of time)    
    umean_realiz = np.zeros(ny, nsavings)

    # Initialize arrays for TTBL thickness parameters
    disp_t = np.zeros(nsavings)   # displacement thickness, delta*
    mom_t  = np.zeros(nsavings)   # momentum     thickness, theta

    """
    Do loop from the first saving of umean (excluding the IC) to the last one, with increment ioutput_cf
    that is read from the input file 'input.i3d'.
    """
    for j in range(ts1, tsn, ioutput_cf):
        
        # Do loop over different realizations
        for i in range(1, nr + 1, 1):
               
            # Read of 'umean' data from 'data/umean' folder
            umean[:,i] = np.loadtxt(f'data_r{i:01d}/umean/umean-ts{j:07d}.txt', skiprows=1, delimiter=None, dtype=np.float64)
            
            # Summing into a sum array for different realizations
            umean_realiz[:,j] = umean_realiz[:,j] + umean[:,i]
    
        #!--- Calculation of thickness parameters ---!

        # Calculate the displacement thickness delta*
        int1 = umean_realiz[:,j]/uwall  # 'integrand 1' 

        # Interpolation at the 6th order of accuracy with a spline of 5th order
        spl = InterpolatedUnivariateSpline(yp, int1, k=5)
        disp_t[j] = spl.integral(y0, yn)

        # Calculate the momentum thickness theta
        int2 = int1 - int1**2  # 'integrand 2' 

        # Interpolation at the 6th order of accuracy with a spline of 5th order
        spl = InterpolatedUnivariateSpline(yp, int2, k=5)
        mom_t[j] = spl.integral(y0, yn)
    
# Save to file

# Plot (optionally) (we can call this function before cf_monitoring so we can read data later and plot with different available quantities)





