#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=wavePacket3D
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097
#SBATCH --output=wavePacket3D.o%j
#SBATCH --error=wavePacket3D.e%j

set -x

# Set number of processors (product must be equal to number of tasks).
ntasks=128
nprocx=64
nprocy=2

# Limit stacksize (adjust to your programs need and core file size).
ulimit -s 204800
ulimit -c 0

# Set OpenMPI configuration.
export OMPI_MCA_osc="ucx"
export OMPI_MCA_pml="ucx"
export OMPI_MCA_btl="self"
export UCX_HANDLE_ERRORS="bt"
export OMPI_MCA_pml_ucx_opal_mem_hooks=1

# Set Intel MPI configuration.
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Define directories.
dirHome=/home/b/b381733/dissertation
dirScratch=/scratch/b/b381733/dissertation
dirNam=${dirHome}/pinc/tests/input
exe=${dirHome}/pinc/bin/pinc
dirWork=${dirScratch}/pinc/tests/wavePacket3D
mkdir ${dirWork}

# Go to work directory.
cd ${dirWork}
rm *

# Copy namelist.
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${dirNam}/input_wavePacket3D.f90 > input.f90

# Run the model.
srun -l --cpu_bind=verbose --hint=nomultithread \
  --distribution=block:cyclic ${exe} 1>run.log 2>&1

exit 0
