#!/bin/bash
#SBATCH --partition=compute
##SBATCH --partition=interactive
#SBATCH --job-name=wkb-wp-1d
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
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
srun --cpu_bind=verbose --distribution=block:cyclic julia --project=examples examples/wavepackets/wkb-wp-1d.jl 1 1>wkb-wp-1d.log 2>&1

# Run the model on interactive partition.
# mpiexec -n 1 julia examples/wavepackets/wkb-wp-1d.jl 1 1>wkb-wp-1d.log 2>&1

exit 0
