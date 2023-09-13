#!/bin/bash
#SBATCH --partition=compute
#SBATCH --account=bb1097
#SBATCH --job-name=barLC
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=49
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097
#SBATCH --output=my_job.%j.out

set -x

#missing ulimit statement
# set Intel MPI environment variables
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Set number of processors (product must be equal to number of tasks).
ntasks=147
nprocx=7
nprocy=21

# Set OpenMP configuration.
export OMP_NUM_THREADS=1

userName=$(whoami)
echo userName

runName=barLC
echo "runName " $runName

inputFile=input_barLC.f90
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
#mpirun -np ${ntasks} ${exe} 1>run.log 2>&1
srun -l --cpu_bind=verbose --hint=nomultithread \
     --distribution=block:cyclic ${exe} 1>run.log 2>&1

exit 0
