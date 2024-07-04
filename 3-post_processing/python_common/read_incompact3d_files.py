#!-----------------------------------------------------------!
#! Here we store functions to read Incompact3d files:        !
#! 'input.i3d', 'post.prm', 'yp.dat'.                        !
#!-----------------------------------------------------------!

import numpy as np

# Define 'read_input_files' function, that reads parameters of simulation 
# from 'input.i3d' and 'post.prm' files.

def read_input_files(filename1,filename2):
    
    # Opening of 'input.i3d' file
    with open(filename1, 'r') as file:
        
        # Read all lines into a list
        lines = file.readlines()
    
        # Extract itype, nx, nz, Lx, Ly, Lz, Re, iswitch_wo 
        # As always, index is 1 less of the line number (Python convention)
        itype      = lines[7]  
        nx         = lines[14]
        nz         = lines[16]
        Lx         = lines[21]
        Ly         = lines[22]
        Lz         = lines[23]
        re         = lines[26]
        iswitch_wo = lines[87]
    
        # Removing characters in front of the extracted strings and the comments:
        # 1) split: the string is split when the specified character is encountered; 
        # 2) we select the portion of string with index inside square brackets;
        # 3) strip: removes leading or trailing whitespaces from the string. 
    
        itype = itype.split('=')[-1].strip()
    
        nx         = nx.split('!')[0]
        nx         = nx.split('=')[-1].strip()
    
        nz         = nz.split('!')[0]
        nz         = nz.split('=')[-1].strip()
    
        Lx         = Lx.split('!')[0]
        Lx         = Lx.split('=')[-1].strip()
    
        Ly         = Ly.split('!')[0]
        Ly         = Ly.split('=')[-1].strip()
    
        Lz         = Lz.split('!')[0]
        Lz         = Lz.split('=')[-1].strip()
    
        re         = re.split('!')[0]
        re         = re.split('=')[-1].strip()
    
        iswitch_wo = iswitch_wo.split('!')[0]
        iswitch_wo = iswitch_wo.split('=')[-1].strip()
    
        # Convert to needed variable type (integer, float, etc.)
        itype      = int(itype)
        nx         = int(nx)
        nz         = int(nz)
        Lx         = np.float64(Lx)
        Ly         = np.float64(Ly)
        Lz         = np.float64(Lz)
        re         = np.float64(re)
        iswitch_wo = int(iswitch_wo)
    
    # Opening of 'post.prm' file
    with open(filename2, 'r') as file:
        
        # Read all lines into a list
        lines = file.readlines()
   
        # Extract needed lines  
        add_string = lines[3]   # Flow case name
        file1      = lines[7]   # First snapshot index
        filen      = lines[8]   # Final snapshot index
        icrfile    = lines[9]   # File increment
        nr         = lines[10]  # Number of flow realizations
        
        # Extract the needed variables
        file1   =   file1.split('#')[0].strip()
        filen   =   filen.split('#')[0].strip()
        icrfile = icrfile.split('#')[0].strip()
        nr      =      nr.split('#')[0].strip()
        
        add_string = add_string.split('!')[0]
        add_string = add_string.rstrip()
        
        # Convert to needed variable type (integer)
        file1      = int(file1)
        filen      = int(filen)
        icrfile    = int(icrfile)
        nr         = int(nr)
        
    # Return to main program with extracted parameters
    return itype, nx, nz, Lx, Ly, Lz, re, iswitch_wo, file1, filen, icrfile, nr, add_string
    

    
    
    
