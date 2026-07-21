#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=mountain_wave
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --hint=nomultithread
#SBATCH --time=0-08:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097
#SBATCH --array=1-6

set -euo pipefail
set -x

export HDF5_USE_FILE_LOCKING=FALSE
export ROMIO_LUSTRE_LOCKING=0
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

case "${SLURM_ARRAY_TASK_ID}" in
  1) UINI=7.0
     RUN=2606_02
     UAMP=3;;
  2) UINI=9.0
     RUN=2606_03
     UAMP=5;;
  3) UINI=11.0
     RUN=2606_03
     UAMP=7;;
  4) UINI=13.0
     RUN=2606_04
     UAMP=5;;
  5) UINI=15.0
     RUN=2606_05
     UAMP=3;;
  6) UINI=17.0
     RUN=2606_06
     UAMP=1;;
  *)
    echo "Invalid SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
    ;;
esac

echo "Running with run ${RUN} with u_ini = ${UINI} and u_amp = ${UAMP}"

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
     julia --project examples/scripts/mountain_wave_ens.jl \
     128 1 1 1 ${UINI} ${RUN} ${UAMP}\
     > mountain_wave_ens_${RUN}.log 2>&1