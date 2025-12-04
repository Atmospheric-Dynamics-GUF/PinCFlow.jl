#!/bin/bash
#SBATCH --partition=general1
#SBATCH --job-name=wkb_wave_packet
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=27
#SBATCH --hint=nomultithread
#SBATCH --time=0-00:15:00
#SBATCH --mail-type=FAIL

set -x

# Configure MPI and HDF5.
julia --project=examples -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
julia --project=examples -e 'using HDF5; HDF5.API.set_libraries!("/home/atmodynamics/public/hdf5-1.14.4-3/src/.libs/libhdf5.so", "/home/atmodynamics/public/hdf5-1.14.4-3/hl/src/.libs/libhdf5_hl.so")'

# Run the model.
mpiexec -n 27 julia examples/scripts/wkb_wave_packet.jl 3 3 3 1>wkb_wave_packet.log 2>&1

exit 0
