#!/bin/bash
#SBATCH --partition=test
#SBATCH --job-name=wavePacket3D
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2600
#SBATCH --mail-type=FAIL
#SBATCH --time=02:00:00

set -x

# cpu settings
ntasks=128
nprocx=64
nprocy=2

# OpenMP settings
export OMP_NUM_THREADS=1

# gprof for MPI
export GMON_OUT_PREFIX=gmon.out-

# env variable export
# export LD_LIBRARY_PATH=LD_LIBRARY_PATH:/home/atmodynamics/jochum/libraries/spack/opt/spack/linux-scientific7-haswell/gcc-4.8.5/hypre-2.26.0-fgfu4ncdxl7ullqryad7jtdzgzqfacjw/lib/

dirHome=/home/atmodynamics/jochum/dissertation
dirScratch=/scratch/atmodynamics/jochum/dissertation

dirNam=${dirHome}/pinc/tests/input
exe=${dirHome}/pinc/bin/pinc
dirWork=${dirScratch}/pinc/tests/wavePacket3D

mkdir ${dirWork}

# go to work

cd ${dirWork}
rm *

# copy namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${dirNam}/input_wavePacket3D.f90 > input.f90

# run the raytracer
mpirun -np ${ntasks} ${exe} 1>run.log 2>&1

exit 0
