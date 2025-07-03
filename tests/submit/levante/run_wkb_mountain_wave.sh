#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=wkb_mountain_wave
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --hint=nomultithread
#SBATCH --time=0-00:10:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set number of processors (product must be equal to number of tasks).
ntasks=64
nprocx=8
nprocy=8

# Define directories.
user=$(whoami)
input=/home/b/${user}/pinc/tests/input
binary=/home/b/${user}/pinc/bin
scratch=/scratch/b/${user}/pinc/tests/wkb_mountain_wave
mkdir -p ${scratch}

# Go to work directory.
cd ${scratch}/ && rm -r *

# Copy namelist.
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${input}/input_wkb_mountain_wave.f90 > input.f90

# Copy the binary.
cp ${binary}/pinc .

# Run the model.
mpiexec -n ${ntasks} ./pinc 1>run.log 2>&1

exit 0
