#!/bin/bash
#SBATCH --job-name=3DWP
#SBATCH --partition=compute
##SBATCH --partition=interactive
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097
#SBATCH --output=my_job.%j.out

#set -x
set -e

# limit stacksize ... adjust to your programs need
# and core file size
ulimit -s 204800
ulimit -c 0

# set Intel MPI environment variables
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Set number of processors (product must be equal to number of tasks).
#ntasks=150
#nprocx=150
#nprocy=1
ntasks=64
nprocx=64
nprocy=1

# Set OpenMP configuration.
export OMP_NUM_THREADS=1

userName=b382140
echo userName

runName=3DWP_LES
echo "runName " $runName

inputFile=input_3DWP_LES.f90
echo "inputFile "$inputFile

# Define directories.
dirHome=/home/b/${userName}/pinc-flow-tracer
dirScratch=/scratch/b/${userName}/output/${runName}
dirSaveCode=/home/b/${userName}/output/${runName}

dirInput=${dirHome}/input
dirSubmit=${dirHome}/submit
exe=${dirHome}/bin/pinc

mkdir -p ${dirScratch}
mkdir -p ${dirSaveCode}

# Copy source dir, input.f90 in dirSaveCode
cp -rp ${dirHome}/src ${dirSaveCode}/.
cp  -p $BASH_SOURCE   ${dirSaveCode}/.

# Go to work directory.
cd ${dirScratch}/ #&& rm *

# Copy namelist.
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${dirInput}/${inputFile} > input.f90

# Copy script in dirSaveCode
if [ ${dirSaveCode} != ${dirScratch} ]; then
   cp -p ${dirScratch}/input.f90 ${dirSaveCode}/.
fi   

# Run the model.
#mpirun -np ${ntasks} ${exe} 1>run.log 2>&1

srun -l --cpu_bind=verbose --hint=nomultithread\
  --distribution=block:cyclic ${exe} 1>run.log 2>&1

exit 0

