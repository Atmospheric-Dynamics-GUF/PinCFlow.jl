#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=wave_packet
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --hint=nomultithread
#SBATCH --time=01:00:00
#SBATCH --account=bb1097
#SBATCH --mem=0
#SBATCH --output=wave_packet.%j.out
#SBATCH --error=wave_packet.%j.err

set -euo pipefail
set -x

JULIA="$HOME/opt/julia-1.10/bin/julia"

module purge
module load intel-oneapi-mpi/2021.5.0-intel-2021.5.0

# Make Intel MPI use Slurm PMI2 on Levante
export I_MPI_PMI=pmi2
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so

# Robust scratch (SCRATCH may be unset)
SCRATCH_DIR="${SCRATCH:-${TMPDIR:-/tmp/$USER/$SLURM_JOB_ID}}"
mkdir -p "$SCRATCH_DIR/julia_depot"
export JULIA_DEPOT_PATH="$HOME/.julia:$SCRATCH_DIR/julia_depot"

# Configure Julia MPI.jl to use the system Intel MPI
"$JULIA" --project=examples -e 'using MPIPreferences; MPIPreferences.use_system_binary(; library_names=["/sw/spack-levante/intel-oneapi-mpi-2021.5.0-mrcss7/mpi/2021.5.0/lib/release/libmpi.so"])'

mkdir -p examples/results

srun --mpi=pmi2 --cpu-bind=cores \
  "$JULIA" --project=examples examples/scripts/research/multiple_wkb_wave_packet.jl 4 4 4
