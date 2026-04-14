#!/bin/bash
#SBATCH --partition=compute
##SBATCH --partition=interactive
#SBATCH --job-name=wp-3d
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=2
#SBATCH --hint=nomultithread
#SBATCH --time=0-08:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097
#SBATCH --output=%x.out

set -x

# Set Intel MPI configuration on compute partition.
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
srun --distribution=block:cyclic julia --project=examples -e 'using PinCFlow, CairoMakie; wp_3d(; npx = 16, npy = 2, npz = 1)' &> wp-3d.log

# Run the model on interactive partition.
# mpiexec -n 27 julia examples/scripts/wp-3d.jl 3 3 3 1>wp-3d.log 2>&1

exit 0
