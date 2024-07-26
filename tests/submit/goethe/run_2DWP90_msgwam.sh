#!/bin/bash
#SBATCH --partition=general2
#SBATCH --job-name=2DWP90
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=2600
#SBATCH --time=02:00:00

set -x

# no. of processors ntasks must be nprocx * nprocy
ntasks=16
nprocx=16
nprocy=1

userName=$(whoami)
echo userName

runName=2DWP90_msgwam
echo "runName " $runName

inputFile=input_2DWP90_msgwam.f90
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

# Run the model.
mpirun -np ${ntasks} ${exe} 1>run.log 2>&1

exit 0