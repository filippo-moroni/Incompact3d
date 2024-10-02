#!/bin/bash
#SBATCH --job-name=test           # job name
#SBATCH --qos=nolimits            # QOS
#SBATCH --account=ecosister       # account name
#SBATCH --partition=ecosister     # partition name
#SBATCH --time=5-23:59:59         # format HH : MM : SS
#SBATCH --nodes=4                 # node number
#SBATCH --ntasks-per-node=52      # tasks out of 52
#SBATCH --cpus-per-task=1         # CPUs per task
#SBATCH --mem=490g                # memory per node
#SBATCH --exclusive               # exclusive use of the nodes
#SBATCH --exclude=cnode[06,17,19] # to exclude some nodes
#------------------- End of SLURM directives -----------------#

#!--- Description ---!
# This is a .sh file for running 4 different flow realizations
# with incompact3d in ARIES.

# Loading of modules (Intel)
module load intel-2023.1/intel_Comps-2023.1
module load intel-2023.1/intel_MPI-2023.1
module load intel-2023.1/intel_MKL-2023.1

# Calculate total number of tasks (verified that it works)
TOTAL_TASKS=$((SLURM_NTASKS_PER_NODE * SLURM_NNODES))

# Launching
mpirun -np $TOTAL_TASKS ../../build/bin/xcompact3d > out_r1

# Renaming /data folder
mv data data_r1

# Launching
mpirun -np $TOTAL_TASKS ../../build/bin/xcompact3d > out_r2

# Renaming /data folder
mv data data_r2

# Launching
mpirun -np $TOTAL_TASKS ../../build/bin/xcompact3d > out_r3

# Renaming /data folder
mv data data_r3

# Launching
mpirun -np $TOTAL_TASKS ../../build/bin/xcompact3d > out_r4

# Renaming /data folder
mv data data_r4

