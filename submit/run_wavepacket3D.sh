#!/bin/bash
#SBATCH --partition=test
#SBATCH --job-name=test
##switch off multi-threading
##BATCH --hint=nomultithread
#SBATCH --ntasks=1
##usually nodes determined by machine. However, if more memory is needed ...
##SBATCH --nodes=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2600
#SBATCH --mail-type=FAIL
#SBATCH --time=01:00:00

set -x

# no. of processors ntasks must be nprocx * nprocy
ntasks=1
nprocx=1
nprocy=1

# OpenMP settings
export OMP_NUM_THREADS=1

# gprof for MPI
export GMON_OUT_PREFIX=gmon.out-

# env variable export
# export LD_LIBRARY_PATH=LD_LIBRARY_PATH:/home/atmodynamics/jochum/spack/opt/spack/linux-scientific7-x86_64/intel-18.0.3/hypre-2.15.1-j7lo2mzfhd7bnrcivaf6bsaqjwimsp3m/lib/

dirHome=/home/irmgard/Documents/Master/WS22/Masterarbeit/pinc-flow-tracer
dirScratch=/home/irmgard/Documents/Master/WS22/Masterarbeit/output-tracer

dirNam=${dirHome}/input
exe=${dirHome}/bin/pinc
dirWork=${dirScratch}/wavePacket3D

mkdir ${dirWork}

# go to work

cd ${dirWork} && rm *

# copy namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${dirNam}/input_wavePacket3D.f90 > input.f90

# run the raytracer
mpirun -np ${ntasks} ${exe} 1>run.log 2>&1

exit 0
