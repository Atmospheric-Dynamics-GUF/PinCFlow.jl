#!/bin/bash
#SBATCH --partition=test
#SBATCH --job-name=wavePacket3D
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2600
#SBATCH --mail-type=FAIL
#SBATCH --time=02:00:00

set -x

# Set number of processors (product must be equal to number of tasks).
ntasks=128
nprocx=64
nprocy=2

# Set OpenMP configuration.
export OMP_NUM_THREADS=1

# Define directories.
dirHome=/home/atmodynamics/jochum/dissertation
dirScratch=/scratch/atmodynamics/jochum/dissertation
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
mpirun -np ${ntasks} ${exe} 1>run.log 2>&1

exit 0
