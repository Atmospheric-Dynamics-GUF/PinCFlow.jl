#!/bin/bash

# This script must be executed with "source".

# Load the cluster modules that provide the MPI and HDF5 backends.
module load intel-oneapi-mpi/2021.5.0-intel-2021.5.0
module load hdf5/1.12.1-intel-oneapi-mpi-2021.5.0-intel-2021.5.0

# Configure MPI and HDF5.
julia --project=examples/PinCFlowExamples.jl -e 'using MPIPreferences; MPIPreferences.use_system_binary(; library_names=["/sw/spack-levante/intel-oneapi-mpi-2021.5.0-mrcss7/mpi/2021.5.0/lib/release/libmpi.so"])'
julia --project=examples/PinCFlowExamples.jl -e 'using HDF5; HDF5.API.set_libraries!("/sw/spack-levante/hdf5-1.12.1-jmeuy3/lib/libhdf5.so", "/sw/spack-levante/hdf5-1.12.1-jmeuy3/lib/libhdf5_hl.so")'

# Precompile.
julia --project=examples/PinCFlowExamples.jl -e 'using Pkg; Pkg.precompile()'
