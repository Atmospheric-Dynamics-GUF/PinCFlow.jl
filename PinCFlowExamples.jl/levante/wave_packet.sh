#!/bin/bash
##SBATCH --partition=compute
#SBATCH --partition=interactive
#SBATCH --job-name=wave_packet
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=27
#SBATCH --hint=nomultithread
#SBATCH --mem=40G
#SBATCH --time=0-01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
# export I_MPI_PMI=pmi
# export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
# srun --cpu_bind=verbose --distribution=block:cyclic julia --project=PinCFlowExamples.jl -e 'using PinCFlowExamples; wave_packet(; npx = 3, npy = 3, npz = 3)' 1>wave_packet.log 2>&1

# Run the model on interactive partition.
mpiexec -n 27 julia --project=PinCFlowExamples.jl -e 'using PinCFlowExamples; wave_packet(; npx = 3, npy = 3, npz = 3)' 1>wave_packet.log 2>&1

exit 0
