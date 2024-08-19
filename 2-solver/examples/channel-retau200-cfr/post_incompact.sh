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

# Loading of modules
module load gcc-12.1.1/gcc
module load gcc-12.1.1/ompi-4.1.5_nccl

# Calculate total number of tasks
TOTAL_TASKS=$((SLURM_NTASKS_PER_NODE * SLURM_NNODES))

# Copying the needed post.prm file
cp post_1.prm post.prm

# Launching for mean statistics
mpirun -np $TOTAL_TASKS ../../../3-post_processing/post_incompact3d/build/bin/post_incompact3d

# Copying the needed post.prm file
cp post_2.prm post.prm

# Launching again for correlations and TKE budgets
mpirun -np $TOTAL_TASKS ../../../3-post_processing/post_incompact3d/build/bin/post_incompact3d

# Deleting post.prm file
rm post.prm




