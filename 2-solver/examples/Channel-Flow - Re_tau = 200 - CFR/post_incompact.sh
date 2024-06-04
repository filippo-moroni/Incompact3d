#!/bin/bash
#SBATCH --account=fluids          # account name
#SBATCH --partition=high          # partition name
#SBATCH --time=23:59:59           # format HH : MM : SS
#SBATCH --nodes=1                 # node number
#SBATCH --ntasks-per-node=52      # tasks out of 52
#SBATCH --mem=400g                # memory per node
#SBATCH --exclusive               # do not share nodes
#SBATCH --job-name=test_Channel   # job name
#------------------- End of SLURM directives -----------------#

# Loading of modules
module load gcc-12.1.1/gcc
module load gcc-12.1.1/ompi-4.1.5_nccl

# Launching
mpirun -np 52 ../../../3-post_processing/post_incompact3d/build/bin/post_incompact3d
