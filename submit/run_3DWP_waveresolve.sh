#!/bin/bash
#SBATCH --partition=general2
#SBATCH --job-name=3DWP_WR
#SBATCH --ntasks=512
#SBATCH --mem-per-cpu=2600
#SBATCH --time=96:00:00

set -x

# no. of processors ntasks must be nprocx * nprocy
ntasks=512
nprocx=256
nprocy=2

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
dirWork=${dirScratch}/output/2024/2024-01-09/3DWP/3DWP_waveresolve_900min_noprime

mkdir ${dirWork}

# go to work

cd ${dirWork} && rm *

# copy namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${dirNam}/input_3DWP_waveresolve.f90 > input.f90

# run the raytracer
mpirun -np ${ntasks} ${exe} 1>run.log 2>&1

exit 0
