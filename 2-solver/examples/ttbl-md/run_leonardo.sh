#!/bin/bash
#SBATCH --job-name=test           # job name
#SBATCH --account= ...            # account name
#SBATCH --partition=dcgp_usr_prod # partition name
#SBATCH --time=23:59:59           # format HH : MM : SS
#SBATCH --nodes=1                 # node number
#SBATCH --ntasks-per-node=112     # tasks out of 112
#SBATCH --cpus-per-task=1         # CPUs per task
#SBATCH --mem=480g                # memory per node
#SBATCH --exclusive               # exclusive use of the nodes
#------------------- End of SLURM directives -----------------#

#!--- Description ---!
# This is a .sh file for running with incompact3d in 
# Leonardo Data Centric General Purpose (DCGP).

# Loading of modules
module load ...

# Calculate total number of tasks
TOTAL_TASKS=$((SLURM_NTASKS_PER_NODE * SLURM_NNODES))

# Launching
mpirun -np $TOTAL_TASKS ../../build/bin/xcompact3d > out
