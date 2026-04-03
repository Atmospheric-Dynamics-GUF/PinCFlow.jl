#!/bin/bash
##SBATCH --partition=compute
#SBATCH --partition=interactive
#SBATCH --job-name=hot_bubble
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
# srun --cpu_bind=verbose --distribution=block:cyclic julia --project=examples/PinCFlowExamples.jl -e 'using PinCFlowExamples; hot_bubble()' 1>hot_bubble.log 2>&1

# Run the model on interactive partition.
mpiexec -n 9 julia --project=examples/PinCFlowExamples.jl -e 'using PinCFlowExamples; hot_bubble()' 1>hot_bubble.log 2>&1

exit 0
