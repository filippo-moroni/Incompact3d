
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



