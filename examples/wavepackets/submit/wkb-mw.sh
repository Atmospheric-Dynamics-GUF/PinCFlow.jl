#!/bin/bash
#SBATCH --partition=compute
##SBATCH --partition=interactive
#SBATCH --job-name=wkb-mw-3d
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
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
srun --distribution=block:cyclic julia --project=examples examples/wavepackets/wkb-mw.jl 4 4 1 1>wkb-mw-3d.log 2>&1

# Run the model on interactive partition.
# mpiexec -n 16 julia --project=examples examples/wavepackets/wkb-mw.jl 4 4 1 1>wkb-mw.log 2>&1

exit 0
