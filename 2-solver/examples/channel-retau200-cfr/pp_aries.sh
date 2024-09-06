#!/bin/bash
#SBATCH --job-name=post           # job name
#SBATCH --account=fluids          # account name
#SBATCH --partition=high          # partition name
#SBATCH --time=23:59:59           # format HH : MM : SS
#SBATCH --nodes=1                 # node number
#SBATCH --ntasks-per-node=52      # tasks out of 52
#SBATCH --cpus-per-task=1         # CPUs per task
#SBATCH --mem=490g                # memory per node
#SBATCH --exclusive               # exclusive use of the nodes
#SBATCH --exclude=cnode[06,17,19] # to exclude some nodes
#------------------- End of SLURM directives -----------------#

#!--- Description ---!
# This is a .sh file for post-processing with post_incompact3d in ARIES.

# Loading of modules (Intel)
module load intel-2023.1/intel_Comps-2023.1
module load intel-2023.1/intel_MPI-2023.1
module load intel-2023.1/intel_MKL-2023.1

# Calculate total number of tasks (verified that it works)
TOTAL_TASKS=$((SLURM_NTASKS_PER_NODE * SLURM_NNODES))

# Copying the needed post.prm file 
cp post_1.prm post.prm

# Launching for mean statistics
mpirun -np $TOTAL_TASKS ../../../3-post_processing/post_incompact3d/build/bin/post_incompact3d > out_pp_mean_stats

# Copying the needed post.prm file 
cp post_2.prm post.prm

# Launching again for correlations and TKE budgets
mpirun -np $TOTAL_TASKS ../../../3-post_processing/post_incompact3d/build/bin/post_incompact3d > out_pp_extra_stats


