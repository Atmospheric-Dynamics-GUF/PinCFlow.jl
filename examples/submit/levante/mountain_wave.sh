#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=mountainwave
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=75
#SBATCH --time=0-00:10:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bb1097

set -x

# Set number of processors (product must be equal to number of tasks).
ntasks=75
nprocx=75
nprocy=1

# Set MPI binary.
mpiexec=`julia -e 'using MPICH_jll; println(MPICH_jll.mpiexec_path)'`

# Define directories.
dirScratch=/scratch/b/b381733/dissertation/pincflow/examples/mountain_wave
dirHome=/home/b/b381733/dissertation/code/pincflow
dirCode=${dirHome}/src
dirSubmit=${dirHome}/examples/submit

# Create work directory and go to it.
mkdir -p ${dirScratch}
cd ${dirScratch}/ && rm -r *

# Copy the code and run script.
cp -r ${dirCode} .
cp ${dirSubmit}/mountain_wave.jl .

# Run the model.
${mpiexec} -n ${ntasks} julia mountain_wave.jl ${nprocx} ${nprocy} 1>run.log 2>&1

exit 0