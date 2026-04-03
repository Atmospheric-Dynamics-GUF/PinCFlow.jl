#!/bin/bash
##SBATCH --partition=compute
#SBATCH --partition=interactive
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
# export I_MPI_PMI=pmi
# export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
# srun --cpu_bind=verbose --distribution=block:cyclic julia --project=examples/PinCFlowExamples.jl -e 'using PinCFlowExamples; wkb_mountain_wave()' 1>wkb_mountain_wave.log 2>&1

# Run the model on interactive partition.
mpiexec -n 27 julia --project=examples/PinCFlowExamples.jl -e 'using PinCFlowExamples; wkb_mountain_wave()' 1>wkb_mountain_wave.log 2>&1

exit 0
