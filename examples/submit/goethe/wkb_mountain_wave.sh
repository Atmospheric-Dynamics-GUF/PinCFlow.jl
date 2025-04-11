#!/bin/bash
#SBATCH --partition=general1
#SBATCH --job-name=wkb_mountain_wave
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --hint=nomultithread
#SBATCH --time=0-00:10:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set the number of tasks.
ntasks=20

# Define the work directory.
user=$(whoami)
scratch=/scratch/atmodynamics/${user}/pinc/examples/wkb_mountain_wave
mkdir -p ${scratch}

# Configure MPI and HDF5.
julia --project=../../../ -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
julia --project=../../../ -e 'using HDF5; HDF5.API.set_libraries!("/home/atmodynamics/public/hdf5-1.14.4-3/src/.libs/libhdf5.so", "/home/atmodynamics/public/hdf5-1.14.4-3/hl/src/.libs/libhdf5_hl.so")'

# Run the model.
mpiexec -n ${ntasks} julia --project=../../../ --check-bounds=no --math-mode=fast wkb_mountain_wave.jl 1>${scratch}/run.log 2>&1

exit 0