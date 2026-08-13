#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=waveBoussinesq
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=8
#SBATCH --hint=nomultithread
#SBATCH --mem=40G
#SBATCH --time=0-01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Run the model on compute partition.
srun --distribution=block:cyclic julia --project=examples -e 'using PinCFlow, CairoMakie; wave_Boussinesq(; 
x_size = 320,
y_size = 32,
z_size = 320,
npx = 32, 
npy = 2,
output_interval = 360.0, 
tmax = 3600.0, 
output_file = "wave_Boussinesq.h5",
)' &> wave_Boussinesq.log

exit 0
