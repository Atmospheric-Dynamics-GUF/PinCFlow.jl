#!/bin/bash
#SBATCH --partition=general1
#SBATCH --job-name=mountain_wave
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --hint=nomultithread
#SBATCH --time=0-00:10:00
#SBATCH --mail-type=FAIL

set -x

# Configure MPI and HDF5.
julia --project=examples -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
julia --project=examples -e 'using HDF5; HDF5.API.set_libraries!("/home/atmodynamics/public/hdf5-1.14.4-3/src/.libs/libhdf5.so", "/home/atmodynamics/public/hdf5-1.14.4-3/hl/src/.libs/libhdf5_hl.so")'

# Run the model.
mpiexec -n 64 julia --project examples/submit/mountain_wave.jl 8 8 1>mountain_wave.log 2>&1

exit 0
