#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=cold_bubble
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --hint=nomultithread
#SBATCH --mem=90G
#SBATCH --time=0-01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
srun --distribution=block:cyclic julia --project=examples -e 'using PinCFlow, CairoMakie; cold_bubble(; npx = 8, npz = 8)' &> cold_bubble.log

exit 0
