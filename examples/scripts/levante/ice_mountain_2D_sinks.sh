#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=mountain_wave
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --hint=nomultithread
#SBATCH --time=0-04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097
#SBATCH --array=1

set -euo pipefail
set -x

export HDF5_USE_FILE_LOCKING=FALSE
export ROMIO_LUSTRE_LOCKING=0
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so



case "${SLURM_ARRAY_TASK_ID}" in
  0) TAUQ=0.0
     TAUQV=0.0
     RUN=2801_;;
  1) TAUQ=0.00000000001
     TAUQV=3000.0
     RUN=2801_02;;
  2) TAUQ=0.0
     TAUQV=3000.0
     RUN=2801_11;;
  3) TAUQ=0.00000000001
     TAUQV=0.0
     RUN=2801_12;;
  *)
    echo "Invalid SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
    ;;
esac

echo "Running with tau_q_sink = ${TAUQ}, tau_qv_source = ${TAUQV}, run = ${RUN}"

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
     julia --project exp/ice_mountain_2D.jl \
     8 1 16 1 ${TAUQ} ${TAUQV} ${RUN} \
     > ice_mountain_2D_${RUN}.log 2>&1