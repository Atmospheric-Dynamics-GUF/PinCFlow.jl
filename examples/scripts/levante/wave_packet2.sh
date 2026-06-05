#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=wave_packet
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --hint=nomultithread
#SBATCH --time=0-00:20:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097
##SBATCH --array=0-1

set -euo pipefail
set -x

export HDF5_USE_FILE_LOCKING=FALSE
export ROMIO_LUSTRE_LOCKING=0
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
     julia --project examples/scripts/wave_packet.jl \
     8 4 1 1\
     > wave_packet.log 2>&1