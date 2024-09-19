"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: In this file, we store useful subroutines for post-processing
!              of Temporal Turbulent Boundary Layer (TTBL) simulations.
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

def calculate_ttbl_delta_99(y,mean_u):

    # Initialize the index
    j = 0
    
    # Calculate the index at which the BL thickness delta99 is and delta_99 itself
    while mean_u[j] > mean_u[0]*0.01: 
        
        # Boundary layer thickness delta_99
        bl_thick = y[j]
        
        # Increment the index
        j = j + 1

    # Rename index for better clarity
    bl_thick_j = j
    
    # Return to main program with bl_thickness and its related j index
    return (bl_thick, bl_thick_j)