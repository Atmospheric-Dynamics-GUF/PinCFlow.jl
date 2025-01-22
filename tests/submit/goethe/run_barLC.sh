#!/bin/bash
#SBATCH --partition=general1
#SBATCH --job-name=barLC
#SBATCH --ntasks=196
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=FAIL
#SBATCH --time=0-00:10:00

set -x

# Set number of processors (product must be equal to number of tasks).
ntasks=196
nprocx=7
nprocy=28

userName=$(whoami)
echo userName

runName=barLC
echo "runName " $runName

inputFile=input_barLC.f90
echo "inputFile "$inputFile

# Define directories.
dirHome=/home/atmodynamics/${userName}/PF/pinc
dirScratch=/scratch/atmodynamics/${userName}/PF/runs/${runName}
dirSaveCode=/home/atmodynamics/${userName}/PF/runs/${runName}

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

# Copy the binary.
cp ${exe} .

# Run the model.
mpirun -np ${ntasks} ./pinc 1>run.log 2>&1

exit 0
