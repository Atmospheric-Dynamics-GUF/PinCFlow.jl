#!/bin/bash
##SBATCH --partition=compute
#SBATCH --partition=interactive
#SBATCH --job-name=wkb_mountain_wave
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --hint=nomultithread
#SBATCH --time=0-00:10:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set Intel MPI configuration on compute partition.
# export I_MPI_PMI=pmi
# export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Define the work directory.
user=$(whoami)
scratch=/scratch/b/${user}/pincflow/examples/wkb_mountain_wave
mkdir -p ${scratch}

# Configure MPI and HDF5.
julia --project=examples -e 'using MPIPreferences; MPIPreferences.use_system_binary(; library_names=["/sw/spack-levante/intel-oneapi-mpi-2021.5.0-mrcss7/mpi/2021.5.0/lib/release/libmpi.so"])'
julia --project=examples -e 'using HDF5; HDF5.API.set_libraries!("/sw/spack-levante/hdf5-1.12.1-jmeuy3/lib/libhdf5.so", "/sw/spack-levante/hdf5-1.12.1-jmeuy3/lib/libhdf5_hl.so")'

# Run the model on compute partition.
# srun --cpu_bind=verbose --distribution=block:cyclic julia --project examples/submit/wkb_mountain_wave.jl ${scratch} 1>${scratch}/run.log 2>&1

# Run the model on interactive partition.
mpiexec -n 64 julia --project examples/submit/wkb_mountain_wave.jl ${scratch} 1>${scratch}/run.log 2>&1

exit 0
