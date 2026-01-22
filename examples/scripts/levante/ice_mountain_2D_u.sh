#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=mountain_wave
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --hint=nomultithread
#SBATCH --time=0-03:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -euo pipefail
set -x

export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Julia environment
julia --project -e 'import Pkg; Pkg.instantiate()'

# MPI
julia --project -e '
using MPIPreferences
MPIPreferences.use_system_binary(
    library_names=["/sw/spack-levante/intel-oneapi-mpi-2021.5.0-mrcss7/mpi/2021.5.0/lib/release/libmpi.so"]
)
'

# HDF5
julia --project -e '
using HDF5
HDF5.API.set_libraries!(
    "/sw/spack-levante/hdf5-1.12.1-jmeuy3/lib/libhdf5.so",
    "/sw/spack-levante/hdf5-1.12.1-jmeuy3/lib/libhdf5_hl.so"
)
'

# Run
srun --cpu_bind=verbose \
     julia --project examples/scripts/ice_mountain_2D_u.jl \
     8 1 16 1 \
     > ice_mountain_2D_u.log 2>&1