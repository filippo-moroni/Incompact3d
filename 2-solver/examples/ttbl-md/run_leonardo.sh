#!/bin/bash
#SBATCH --job-name=test           # job name
#SBATCH --account=CNHPC_1572529   # account name
#SBATCH --partition=dcgp_qos_dbg  # partition name
#SBATCH --time=00:29:59           # format HH : MM : SS
#SBATCH --nodes=1                 # node number
#SBATCH --ntasks-per-node=112     # tasks out of 112
#SBATCH --cpus-per-task=1         # CPUs per task
#SBATCH --mem=480g                # memory per node
#SBATCH --exclusive               # exclusive use of the nodes
#------------------- End of SLURM directives -----------------#

#!--- Description ---!
# This is a .sh file for running with incompact3d in 
# Leonardo Data Centric General Purpose (DCGP).

# Loading of modules (Intel)
module load intel-oneapi-compilers/2023.2.1
module load intel-oneapi-mpi/2021.10.0
module load intel-oneapi-mkl/2023.2.0 

# Calculate total number of tasks (verified that it works)
TOTAL_TASKS=$((SLURM_NTASKS_PER_NODE * SLURM_NNODES))

# Launching
mpirun -np $TOTAL_TASKS ../../build/bin/xcompact3d > out
