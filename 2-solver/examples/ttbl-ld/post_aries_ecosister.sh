#!/bin/bash
#SBATCH --job-name=post           # job name
#SBATCH --qos=nolimits            # QOS
#SBATCH --account=ecosister       # account name
#SBATCH --partition=ecosister     # partition name
#SBATCH --time=5-23:59:59         # format DD : HH : MM : SS (max 6 days)
#SBATCH --nodes=1                 # node number
#SBATCH --ntasks-per-node=52      # tasks out of 52
#SBATCH --cpus-per-task=1         # CPUs per task
#SBATCH --mem=490g                # memory per node
#SBATCH --exclusive               # exclusive use of the nodes
#SBATCH --exclude=cnode[06,17,19] # to exclude some nodes
#------------------- End of SLURM directives -----------------#

#!--- Loading of modules (ARIES) ---!

# Intel
#module load intel-2023.1/intel_Comps-2023.1
#module load intel-2023.1/intel_MPI-2023.1
#module load intel-2023.1/intel_MKL-2023.1

# Loading of modules (GNU)
module load gcc-12.1.1/gcc
module load gcc-12.1.1/ompi-4.1.5_nccl

# Calculate total number of tasks
TOTAL_TASKS=$((SLURM_NTASKS_PER_NODE * SLURM_NNODES))

# Copying the needed post.prm file
cp post_1.prm post.prm

# Launching for mean statistics
mpirun -np $TOTAL_TASKS ../../../3-post_processing/post_incompact3d/build/bin/post_incompact3d > out_mean_stats

# Copying the needed post.prm file
cp post_2.prm post.prm

# Launching again for correlations and TKE budgets
mpirun -np $TOTAL_TASKS ../../../3-post_processing/post_incompact3d/build/bin/post_incompact3d > out_extra_stats

# Deleting post.prm file
rm post.prm
