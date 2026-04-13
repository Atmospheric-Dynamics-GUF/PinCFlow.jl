#!/bin/bash
#SBATCH --partition=compute
##SBATCH --partition=interactive
#SBATCH --job-name=mountain_wave
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --hint=nomultithread
#SBATCH --time=0-04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

RUN="0213_01"

# Set Intel MPI configuration on compute partition.
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

project_dir="/home/b/b383844/PinCFlow/PinCFlow/."

julia --project=${project_dir} -e 'import Pkg; Pkg.instantiate()'

# Configure MPI and HDF5.
julia --project=${project_dir} -e 'using MPIPreferences; MPIPreferences.use_system_binary(; library_names=["/sw/spack-levante/intel-oneapi-mpi-2021.5.0-mrcss7/mpi/2021.5.0/lib/release/libmpi.so"])'
julia --project=${project_dir} -e 'using HDF5; HDF5.API.set_libraries!("/sw/spack-levante/hdf5-1.12.1-jmeuy3/lib/libhdf5.so", "/sw/spack-levante/hdf5-1.12.1-jmeuy3/lib/libhdf5_hl.so")' 
# Run the model on compute partition.
srun --cpu_bind=verbose \
     julia --project="$1" "$2" \
     8 1 16 1\
     > ice_mountain_2D_${RUN}.log 2>&1

# Run the model on interactive partition.
#mpiexec -n 4 julia --project=${project_dir} ${1} 4 1 1 1>ice_dump.log 2>&1

exit 0
