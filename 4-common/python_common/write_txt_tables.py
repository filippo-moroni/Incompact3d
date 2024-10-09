"""
!-----------------------------------------------------------------------------!
! DESCRIPTION: Python function to create and write tables.
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
"""

# Libraries
from tabulate import tabulate

def write_txt_tables(f, data_arrays, titles):

    # Loop to create and write multiple tables to the file 'f'
    for title, data in zip(titles, data_arrays):
        f.write(f"!----- {title} -----!\n\n")
        f.write(tabulate(data, headers="firstrow", tablefmt="fancy_grid") + "\n")

# Inside your main code, where you want to create the tables
data_arrays = [data1, data2, data3]
titles = ["Inputs", "Numerics-related parameters", "Outputs"]

with open("sim_settings.txt", "w") as f:
    f.write("!----- Simulation Settings -----!\n\n")
    create_and_write_tables(f, data_arrays, titles)

