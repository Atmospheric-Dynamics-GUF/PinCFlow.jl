#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=wkb_wave_packet
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --hint=nomultithread
#SBATCH --mem=100G
#SBATCH --time=0-01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
srun --distribution=block:cyclic julia --project=examples -e 'using PinCFlow, CairoMakie; wkb_wave_packet(; npx = 4, npy = 4, npz = 4)' &> wkb_wave_packet.log

exit 0
