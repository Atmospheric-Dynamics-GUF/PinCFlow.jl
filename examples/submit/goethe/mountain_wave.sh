#!/bin/bash
#SBATCH --partition=general1
#SBATCH --job-name=mountain_wave
#SBATCH --ntasks=75
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=FAIL
#SBATCH --time=0-00:10:00

set -x

# Set number of processors (product must be equal to number of tasks).
ntasks=75
nprocx=75
nprocy=1

# Define directories.
dirScratch=/scratch/atmodynamics/jochum/dissertation/pincflow/examples/mountain_wave
dirHome=/home/atmodynamics/jochum/dissertation/pincflow
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