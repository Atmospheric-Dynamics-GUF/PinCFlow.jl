#!/bin/bash
#SBATCH --job-name=IGW
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097
#SBATCH --output=my_job.%j.out

# limit stacksize ... adjust to your programs need
# and core file size
ulimit -s 204800
ulimit -c 0

set -x
#set -e

# Set number of processors (product must be equal to number of tasks).
ntasks=1
nprocx=1
nprocy=1

# Set OpenMP configuration.
export OMP_NUM_THREADS=1

userName=$(whoami)
echo userName

runName=IGW
echo "runName " $runName

inputFile=input_IGW.f90
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

