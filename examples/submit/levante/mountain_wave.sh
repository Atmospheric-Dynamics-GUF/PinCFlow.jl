#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=mountain_wave
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --hint=nomultithread
#SBATCH --time=0-00:10:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set number of processors (product must be equal to number of tasks).
ntasks=64
nprocx=8
nprocy=8

# Define directories.
dirScratch=/scratch/b/b381733/dissertation/pinc/examples/mountain_wave
dirHome=/home/b/b381733/dissertation/code/pinc
dirCode=${dirHome}/src
dirSubmit=${dirHome}/examples/submit

# Create work directory and go to it.
mkdir -p ${dirScratch}
cd ${dirScratch}/ && rm -r *

# Copy the code and run script.
cp -r ${dirCode} .
cp ${dirSubmit}/mountain_wave.jl .

# Configure MPI and run the model (system binary).
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary(; library_names=["/sw/spack-levante/intel-oneapi-mpi-2021.5.0-mrcss7/mpi/2021.5.0/lib/release/libmpi.so"])'
mpiexec -n ${ntasks} julia --project=. --check-bounds=no --math-mode=fast mountain_wave.jl ${nprocx} ${nprocy} 1>run.log 2>&1

exit 0