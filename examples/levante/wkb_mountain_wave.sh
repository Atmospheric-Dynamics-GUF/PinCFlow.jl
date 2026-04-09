#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=wkb_mountain_wave
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=27
#SBATCH --hint=nomultithread
#SBATCH --mem=80G
#SBATCH --time=0-01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
srun --distribution=block:cyclic julia --project=examples -e 'using PinCFlow, CairoMakie; wkb_mountain_wave(; npx = 3, npy = 3, npz = 3)' &> wkb_mountain_wave.log

exit 0
