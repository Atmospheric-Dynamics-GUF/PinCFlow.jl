#!/bin/bash
##SBATCH --partition=compute
#SBATCH --partition=interactive
#SBATCH --job-name=periodic_hill
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9
#SBATCH --hint=nomultithread
#SBATCH --mem=15G
#SBATCH --time=0-00:15:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
# export I_MPI_PMI=pmi
# export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
# srun --cpu_bind=verbose --distribution=block:cyclic julia examples/scripts/periodic_hill.jl 4 4 1>periodic_hill.log 2>&1

# Run the model on interactive partition.
mpiexec -n 9 julia examples/scripts/periodic_hill.jl 3 3 1>periodic_hill.log 2>&1

exit 0
