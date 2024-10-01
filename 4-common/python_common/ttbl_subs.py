"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: In this file, we store useful subroutines for post-processing
!              of Temporal Turbulent Boundary Layer (TTBL) simulations.
!              List of subroutines:
!               - extract_re_tau_value;
!               - calculate_ttbl_delta_99.
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this script, we read the friction Reynolds number Re_tau 
!              from .xdmf files of snapshots of TTBL simulations.    
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
from lxml import etree

def extract_re_tau_value(file_path):
    
    # Parse the XDMF file (read)
    tree = etree.parse(file_path)
    root = tree.getroot()
    
    # Find all comments
    comments = tree.xpath('//comment()')
    
    for comment in comments:
        if "Re_tau" in comment.text:
            # Extract the numerical value after "Re_tau ="
            try:
                # Assuming the comment format is "Friction Reynolds number, Re_tau = value"
                re_tau_part = comment.text.split("Re_tau =")[1]
                re_tau = float(re_tau_part.strip().split()[0])
                
                # Round-off Re_tau value
                re_tau = round(re_tau)
                
                # Return Re_tau value to main program
                return re_tau
            except (IndexError, ValueError):
                print("Error extracting Re_tau value from comment.")
                return None
    
    # Print error if Re_tau value is not found
    print("Re_tau comment not found.")
    return None

#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: With this script, we calculate the boundary layer thickness
!              delta_99 for a TTBL with translating wall.
!              See for reference Kozul et al. (2016).    
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
from scipy.interpolate import InterpolatedUnivariateSpline

def calculate_ttbl_thick_params(mean_u,y,uwall):

    """
    !-----------------------------------------------------------!
     TTBL thickness delta_99 calculation (Kozul et al. (2016)). 
    !-----------------------------------------------------------!    
    """
    
    # Initialize the TTBL delta_99 index
    j = 0
    
    # Calculate the index at which the BL thickness delta99 is and delta_99 itself
    while mean_u[j] > mean_u[0]*0.01: 
        
        # Boundary layer thickness delta_99
        delta_99 = y[j]
        
        # Increment the index
        j = j + 1

    # Rename index for better clarity
    delta_99_j = j
    
    """
    !-----------------------------------------------------------!
     TTBL displacement thickness (disp_t) delta* calculation 
     (Kozul et al. (2016)). 
    !-----------------------------------------------------------!    
    """
    
    # First element of 'yp.dat' (y = 0)
    y0 = y[0]

    # Last element of 'yp.dat' (y = Ly, height of the domain)    
    yn = y[-1]
    
    # Calculate the displacement thickness delta* (int1: 'integrand 1')
    int1 = mean_u/uwall 

    # Interpolation at the 6th order of accuracy with a spline of 5th order, 
    # that passes through all data points
    spl = InterpolatedUnivariateSpline(y, int1, k=5)
    disp_t = spl.integral(y0, yn)
    
    """
    !-----------------------------------------------------------!
     TTBL momentum thickness (mom_t) theta calculation 
     (Kozul et al. (2016)). This is the same as a standard
     spatial TBL. 
    !-----------------------------------------------------------!    
    """

    # Calculate the momentum thickness theta (int2: 'integrand 2') 
    int2 = int1 - int1**2 

    # Interpolation at the 6th order of accuracy with a spline of 5th order,
    # that passes through all data points
    spl = InterpolatedUnivariateSpline(y, int2, k=5)
    mom_t = spl.integral(y0, yn)
        
     
    # Return to main program with TTBL thicknesses and delta_99 related j index
    return (delta_99, delta_99_j, disp_t, mom_t)

#!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!


