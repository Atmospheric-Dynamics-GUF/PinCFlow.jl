#!/bin/bash
#SBATCH --partition=compute
##SBATCH --partition=interactive
#SBATCH --job-name=wp-1d
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=2
#SBATCH --hint=nomultithread
#SBATCH --time=0-08:00:00
#SBATCH --mail-type=ALL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
srun --cpu_bind=verbose --distribution=block:cyclic julia examples/scripts/wp-1d.jl 4 4 1 1>wp-1d.log 2>&1

# Run the model on interactive partition.
# mpiexec -n 4 julia examples/scripts/wp-1d.jl 2 2 1 1>wp-1d.log 2>&1

exit 0
