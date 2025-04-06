#!/bin/bash
#SBATCH --partition=general1
#SBATCH --job-name=mountain_wave
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --hint=nomultithread
#SBATCH --time=0-00:10:00
#SBATCH --mail-type=FAIL

set -x

# Set number of processors (product must be equal to number of tasks).
ntasks=20
nprocx=5
nprocy=4

# Define directories.
dirScratch=/scratch/atmodynamics/jochum/dissertation/pinc/examples/mountain_wave
dirHome=/home/atmodynamics/jochum/dissertation/pinc
dirCode=${dirHome}/src
dirSubmit=${dirHome}/examples/submit

# Create work directory and go to it.
mkdir -p ${dirScratch}
cd ${dirScratch}/ && rm -r *

# Copy the code and run script.
cp -r ${dirCode} .
cp ${dirSubmit}/mountain_wave.jl .

# Configure MPI and run the model (system binary).
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
mpiexec -n ${ntasks} julia mountain_wave.jl ${nprocx} ${nprocy} 1>run.log 2>&1

exit 0