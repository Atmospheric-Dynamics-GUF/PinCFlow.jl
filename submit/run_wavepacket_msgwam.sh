#!/bin/bash
#SBATCH --partition=test
#SBATCH --job-name=wavepacket_msgwam
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5200
#SBATCH --mail-type=ALL
#SBATCH --mail-use=s9467794@stud.uni-frankfurt.de
#SBATCH --time=02:00:00

set -x

# no. of processors ntasks must be nprocx * nprocy
ntasks=64
nprocx=64
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
dirWork=${dirScratch}/output/wavepacket_msgwam_smooth

mkdir ${dirWork}

# go to work

cd ${dirWork} && rm *

# copy namelist
sed -e "s/{nprocx}/${nprocx}/" \
    -e "s/{nprocy}/${nprocy}/" \
        ${dirNam}/input_wavepacket_msgwam.f90 > input.f90

# run the raytracer
mpirun -np ${ntasks} ${exe} 1>run.log 2>&1

exit 0
