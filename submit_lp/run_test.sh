#!/bin/bash
#SBATCH --partition=test
#SBATCH --job-name=wavePacket3D
#SBATCH --nodes=2
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2300
#SBATCH --mail-type=FAIL
#SBATCH --time=02:00:00

set -x

# Set number of processors (product must be equal to number of tasks).
ntasks=4 # 
nprocx=4
nprocy=1

# Set OpenMP configuration.
export OMP_NUM_THREADS=1

userName=$(whoami)
echo userName

runName=$1
echo "runName " $runName

if [ -z "$2" ]
  then
      echo "No second argument supplied"
      inputName=$1
else
     inputName=$2
fi

# define namelist file
inputFile=input_${inputName}.f90
#input_3DWP_LES_ice.f90 #does not work
echo "inputFile "$inputFile

# Define directories.
dirHome=/home/${userName}/PF/pinc
dirScratch=/home/${userName}/PF/runs/${runName}
dirSaveCode=/home/${userName}/PF/code_runs/${runName}

dirInput=${dirHome}/input
dirSubmit=${dirHome}/submit
exe=${dirHome}/bin/pinc

mkdir -p ${dirScratch}
mkdir -p ${dirSaveCode}

# save hash latest commit
git rev-parse HEAD > ${dirHome}/src/git_hash

# Copy source dir, input.f90 in dirSaveCode
cp -rp ${dirHome}/src ${dirSaveCode}/.
cp  -p $BASH_SOURCE   ${dirSaveCode}/.

# Go to work directory.
cd ${dirScratch}/ && rm -r *
cp -p ${dirHome}/src/git_hash .

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
