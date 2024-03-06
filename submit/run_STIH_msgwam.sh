#!/bin/bash
#SBATCH --partition=test
#SBATCH --job-name=stihWKB
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=02:00:00

set -x

# no. of processors ntasks must be nprocx * nprocy
ntasks=1
nprocx=1
nprocy=1

# OpenMP settings
export OMP_NUM_THREADS=1

# gprof for MPI
export GMON_OUT_PREFIX=gmon.out-

dirHome=/home/atmodynamics/knop/pinc-flow-tracer
dirScratch=/scratch/atmodynamics/knop

dirNam=${dirHome}/input
exe=${dirHome}/bin/pinc
dirWork=${dirScratch}/output/2024/2024-03-01/STIH_msgwam

mkdir ${dirWork}

# go to work

cd ${dirWork} && rm *

# copy namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${dirNam}/input_STIH_msgwam.f90 > input.f90

# run the raytracer
mpirun -np ${ntasks} ${exe} 1>run.log 2>&1

exit 0
