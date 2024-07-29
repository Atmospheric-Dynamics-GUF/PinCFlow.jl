#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=2DWP90_msgwam
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097
#SBATCH --output=2DWP90_msgwam.o%j
#SBATCH --error=2DWP90_msgwam.e%j

set -x

# Set number of processors (product must be equal to number of tasks).
ntasks=16
nprocx=16
nprocy=1

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

userName=$(whoami)
echo userName

runName=2DWP90_msgwam
echo "runName " $runName

inputFile=input_2DWP90_msgwam.f90
echo "inputFile "$inputFile

# Define directories.
dirHome=/home/b/${userName}/PF/pinc
dirScratch=/scratch/b/${userName}/PF/runs/${runName}
dirSaveCode=/home/b/${userName}/PF/runs/${runName}

dirInput=${dirHome}/input
dirSubmit=${dirHome}/submit
exe=${dirHome}/bin/pinc

mkdir -p ${dirScratch}
mkdir -p ${dirSaveCode}

# Copy source dir, input.f90 in dirSaveCode
cp -rp ${dirHome}/src ${dirSaveCode}/.
cp  -p $BASH_SOURCE   ${dirSaveCode}/.

# Go to work directory.
cd ${dirScratch}/ && rm -r *

# Copy namelist.
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${dirInput}/${inputFile} > input.f90

# Copy script in dirSaveCode
if [ ${dirSaveCode} != ${dirScratch} ]; then
   cp -p ${dirScratch}/input.f90 ${dirSaveCode}/.
fi

# Run the model.
srun -l --cpu_bind=verbose --hint=nomultithread \
  --distribution=block:cyclic ${exe} 1>run.log 2>&1

exit 0
