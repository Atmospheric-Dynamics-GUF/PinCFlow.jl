#!/bin/bash
##SBATCH --partition=compute
#SBATCH --partition=interactive
#SBATCH --job-name=wkb_mountain_wave
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=27
#SBATCH --hint=nomultithread
#SBATCH --time=0-00:15:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
# export I_MPI_PMI=pmi
# export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
# srun --cpu_bind=verbose --distribution=block:cyclic julia examples/scripts/wkb_mountain_wave.jl 4 4 4 1>wkb_mountain_wave.log 2>&1

# Run the model on interactive partition.
mpiexec -n 27 julia examples/scripts/wkb_mountain_wave.jl 3 3 3 1>wkb_mountain_wave.log 2>&1

exit 0
