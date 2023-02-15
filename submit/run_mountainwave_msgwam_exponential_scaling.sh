#!/bin/bash
#SBATCH --partition=general2
#SBATCH --job-name=mountainwave_msgwam_exponential_scaling
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2600
#SBATCH --mail-type=FAIL
#SBATCH --time=03:00:00

set -x

# no. of processors ntasks must be nprocx * nprocy
ntasks=32
nprocx=32
nprocy=1

# OpenMP settings
export OMP_NUM_THREADS=1

# gprof for MPI
# export GMON_OUT_PREFIX=gmon.out-

# env variable export
# export LD_LIBRARY_PATH=LD_LIBRARY_PATH:/home/atmodynamics/jochum/libraries/spack/opt/spack/linux-scientific7-haswell/gcc-4.8.5/hypre-2.26.0-fgfu4ncdxl7ullqryad7jtdzgzqfacjw/lib/

dirHome=/home/atmodynamics/jochum/dissertation
dirScratch=/scratch/atmodynamics/jochum/dissertation

dirNam=${dirHome}/pinc/input
exe=${dirHome}/pinc/bin/pinc
dirWork=${dirScratch}/pinc/mountainwave_msgwam_exponential_scaling

mkdir ${dirWork}

# go to work

cd ${dirWork}
rm *

# copy namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${dirNam}/input_mountainwave_msgwam_exponential_scaling.f90 > input.f90

# run the raytracer
mpirun -np ${ntasks} ${exe} 1>run.log 2>&1

exit 0
