#!/bin/bash
#SBATCH --partition=test
#SBATCH --job-name=3DWP_MS
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=5200
#SBATCH --time=02:00:00

set -x

# no. of processors ntasks must be nprocx * nprocy
ntasks=8
nprocx=8
nprocy=1

# OpenMP settings
export OMP_NUM_THREADS=1

# gprof for MPI
export GMON_OUT_PREFIX=gmon.out-

# env variable export
# export LD_LIBRARY_PATH=LD_LIBRARY_PATH:/home/atmodynamics/jochum/spack/opt/spack/linux-scientific7-x86_64/intel-18.0.3/hypre-2.15.1-j7lo2mzfhd7bnrcivaf6bsaqjwimsp3m/lib/

dirHome=/home/atmodynamics/knop/pinc-flow-tracer
dirScratch=/scratch/atmodynamics/knop

dirNam=${dirHome}/input
exe=${dirHome}/bin/pinc
dirWork=${dirScratch}/output/2023-09-14/gauss/3DWP_msgwam_150min_gauss_noforce

mkdir ${dirWork}

# go to work

cd ${dirWork} && rm *

# copy namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${dirNam}/input_3DWP_msgwam.f90 > input.f90

# run the raytracer
mpirun -np ${ntasks} ${exe} 1>run.log 2>&1

exit 0
